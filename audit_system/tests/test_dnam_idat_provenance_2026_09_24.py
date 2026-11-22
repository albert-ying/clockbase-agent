"""Deterministic miniature SOFTWARE FIXTURES, not empirical IDAT measurements.

The byte generator and tables below exercise integrity contracts only. No R,
normalization, biological simulation, clock, download or model call is made.
"""
import copy
import gzip
import hashlib
import json
from pathlib import Path
import struct
import tempfile
import unittest

import numpy as np
import pandas as pd

from audit_system.report_2026_09_10 import (
    _idat_binary_readback, _idat_literal_pairs, _idat_verification_profile, verify_dnam_idat_provenance,
)


class IdatFixture:
    """A disclosed software fixture: 2 GSMs, 3 probes, 10 addresses, 4 channels."""
    ids = ['GSM100001', 'GSM100002']
    probes = ['cg00000001', 'cg00000002', 'cg00000003']

    def __init__(self, root):
        self.root = Path(root).resolve(); self.out = self.root/'processing'
        self.out.mkdir(parents=True); self.profile = _idat_verification_profile()
        self.files = {}; self.channel_data = {}
        (self.root/'SOFTWARE_FIXTURE.txt').write_text('Disclosed deterministic software fixture; no empirical measurements or minfi execution.\n')
        series = self.root/'series.soft'; samples = self.root/'samples.soft'
        series.write_text('^SERIES = GSE999001\n'+''.join('!Series_sample_id = '+gsm+'\n' for gsm in self.ids))
        pairs, lines = [], []
        for number,gsm in enumerate(self.ids):
            pair = {'sample_id':gsm,'channels':{}}
            lines.extend(['^SAMPLE = '+gsm, '!Sample_platform_id = GPL23976'])
            for channel in ('Grn','Red'):
                filename = gsm+'_ConditionLabel_'+str(number+1)+'_'+channel+'.idat.gz'
                url = 'https://ftp.ncbi.nlm.nih.gov/geo/samples/'+gsm[:-3]+'nnn/'+gsm+'/suppl/'+filename
                pair['channels'][channel] = {'url':url,'https_url':url,'filename':filename,'basename':filename.rsplit('_',1)[0]}
                lines.append('!Sample_supplementary_file = '+url)
            pairs.append(pair)
        samples.write_text('\n'.join(lines)+'\n')
        self.add('series',series); self.add('samples',samples)
        self.plan = {'schema':'clockbase-dnam-idat-plan-v1','study_id':'GSE999001','created_unix':1,
            'profile':self.profile,'pairs':pairs,'source_files':{key:copy.deepcopy(self.files[key]) for key in ('series','samples')},
            'limits':{'compressed_bytes_per_file':128_000_000,'decompressed_bytes_per_file':256_000_000,
                      'compressed_bytes_total':2_000_000_000,'processing_seconds':3600,'samples':64},
            'new_analysis_profile':True,'model_calls':0,'biological_units_inferred':False,'scientific_status':'CANNOT VERIFY'}
        self.plan_path = self.root/'protocol.json'; self.json(self.plan_path,self.plan); self.add('protocol',self.plan_path)
        acquisition_pairs = []
        for number,pair in enumerate(pairs):
            gsm = pair['sample_id']; joined = {'sample_id':gsm,'channels':{}}
            for channel in ('Grn','Red'):
                suffix = gsm+':'+channel; declaration = pair['channels'][channel]
                data = {'address':np.arange(1,11,dtype=np.int32),
                    'Mean':np.arange(1,11,dtype=np.uint16)*10+100*number+(3 if channel=='Red' else 0),
                    'SD':np.arange(1,11,dtype=np.uint16),'NBeads':np.arange(10,20,dtype=np.uint8)}
                header = {'version':3,'nSNPsRead':10,'RedGreen':0,'Barcode':'111111111111',
                          'ChipType':'BeadChip 8x5','position':'R0'+str(number+1)+'C01'}
                self.channel_data[suffix] = {'data':data,'header':header}
                compressed = self.root/'acquisition'/gsm/channel/'attempt-001'/declaration['filename']
                raw = self.root/'acquisition/uncompressed'/declaration['filename'][:-3]
                compressed.parent.mkdir(parents=True); raw.parent.mkdir(parents=True,exist_ok=True)
                binary = self.binary(header,data); raw.write_bytes(binary); compressed.write_bytes(gzip.compress(binary,mtime=0))
                receipt, marker = Path(str(compressed)+'.source.json'), Path(str(raw)+'.json')
                self.json(receipt,{'url':declaration['https_url'],'retrieved_unix':2,
                    'sha256':self.sha(compressed),'size_bytes':compressed.stat().st_size})
                self.json(marker,{'compressed_sha256':self.sha(compressed),'sha256':self.sha(raw),'size_bytes':raw.stat().st_size})
                for kind,path in [('compressed',compressed),('raw',raw),('source_receipt',receipt),('decompression',marker)]:
                    self.add(kind+':'+suffix,path)
                joined['channels'][channel] = {'path':str(raw),'sha256':self.sha(raw),**declaration}
            acquisition_pairs.append(joined)
        self.acquisition = {'schema':'clockbase-dnam-idat-acquisition-v1','study_id':'GSE999001',
            'plan_path':str(self.plan_path),'plan_sha256':self.sha(self.plan_path),'profile':self.profile,
            'files':copy.deepcopy({k:v for k,v in self.files.items() if ':' in k}), 'pairs':acquisition_pairs,
            'compressed_bytes':sum(v['size_bytes'] for k,v in self.files.items() if k.startswith('compressed:')),
            'recorded_unix':3,'new_clock_runs':0,'model_calls':0,'scientific_status':'CANNOT VERIFY'}
        self.acquisition_path = self.root/'acquisition/acquisition.json'
        self.json(self.acquisition_path,self.acquisition); self.add('acquisition',self.acquisition_path)
        self.outputs = {}
        def output(role,name):
            self.outputs[role] = name
            return self.out/name
        for role,name in [('implementation:rworker','worker.R'),('implementation:inputs','inputs.py')]:
            path = self.out/'implementation'/name; path.parent.mkdir(exist_ok=True)
            path.write_text('# Disclosed software fixture; not executed.\n'); self.add(role,path)
        for suffix,record in self.channel_data.items():
            gsm,channel = suffix.split(':')
            self.json(output('channel_header:'+suffix,gsm+'-'+channel+'-header.json'), {'sample_id':gsm,'channel':channel,**record['header']})
            frame = pd.DataFrame({key:record['data'][key] for key in ('Mean','SD','NBeads')},index=record['data']['address'])
            frame.index.name='address'
            self.frame(output('channel_decoded:'+suffix,gsm+'-'+channel+'-decoded.csv.gz'),frame,index=True)
        self.frame(output('manifest_I','manifest-I.csv'),pd.DataFrame([
            [self.probes[0],1,2,'Grn'],[self.probes[1],3,4,'Red']],columns=['probe_id','AddressA','AddressB','Color']))
        self.frame(output('manifest_II','manifest-II.csv'),pd.DataFrame([[self.probes[2],5]],columns=['probe_id','AddressA']))
        self.frame(output('manifest_controls','manifest-controls.csv'),pd.DataFrame([
            [address,kind,'Both','disclosed fixture'] for address,kind in zip(range(6,11),['NEGATIVE','NORM_A','NORM_C','NORM_G','NORM_T'])],
            columns=['Address','Type','Color','ExtendedType']))
        raw_m = np.array([[20,120],[43,143],[50,150]],dtype=float)
        raw_u = np.array([[10,110],[33,133],[53,153]],dtype=float)
        m = np.array([[200,300],[-5,100],[np.nan,400]],dtype=float)
        u = np.array([[100,-10],[150,50],[100,350]],dtype=float)
        beta = np.maximum(m,0)/(np.maximum(m,0)+np.maximum(u,0)+100)
        det = np.array([[0.01,0.010001],[np.nan,0.005],[0.001,np.inf]])
        flags = np.isfinite(det)&(det>0.01); det_na=~np.isfinite(det); beta_na=~np.isfinite(beta)
        mask=flags|det_na|beta_na; matrix=beta.copy();matrix[mask]=np.nan
        self.matrices = {'raw_meth':raw_m,'raw_unmeth':raw_u,'meth_corrected':m,'unmeth_corrected':u,
                         'beta_unmasked':beta,'detection_p':det,'matrix':matrix}
        for role,values in self.matrices.items():
            filename = {'matrix':'beta-qc.csv.gz','meth_corrected':'meth-corrected.csv.gz','unmeth_corrected':'unmeth-corrected.csv.gz',
                        'raw_meth':'raw-meth.csv.gz','raw_unmeth':'raw-unmeth.csv.gz','beta_unmasked':'beta-unmasked.csv.gz','detection_p':'detection-p.csv.gz'}[role]
            frame = pd.DataFrame(values,index=self.probes,columns=self.ids);frame.index.name='probe_id'
            self.frame(output(role,filename),frame,index=True)
        qc=pd.DataFrame({'sample_id':self.ids,'probe_count':3,'detection_failed_count':flags.sum(axis=0),
            'detection_nonfinite_count':det_na.sum(axis=0),'beta_nonfinite_count':beta_na.sum(axis=0),
            'masked_count':mask.sum(axis=0),'observed_count':(~mask).sum(axis=0),
            'raw_median_M':np.median(raw_m,axis=0),'raw_median_U':np.median(raw_u,axis=0)})
        self.frame(output('sample_qc','sample-qc.csv'),qc)
        for role,name in [('session_info','session-info.txt'),('function_source','function-source.R.txt')]:
            output(role,name).write_text('Disclosed software fixture; no R process executed.\n')
        self.processing={'schema':'clockbase-dnam-idat-processing-v1','status':'PASS','profile':self.profile,
            'sample_ids':self.ids,'array':self.profile['array'],'annotation':{'array':self.profile['array'],'annotation':'fixture'},
            'outputs':self.outputs,'compatibility':{'channel_headers':'PASS','address_sets':'PASS','address_count':10,
                'importer_values':'PASS','manifest_mapping':'PASS','required_manifest_addresses':10,'missing_manifest_addresses':0,
                'probe_count':3,'control_count':5,'typeI_count':2,'typeII_count':1,
                'raw_array_ids':['111111111111:R01C01','111111111111:R02C01']},
            'sample_vs_joint_max_abs_delta':dict.fromkeys(self.ids,0),
            'packages':{'minfi':'1.56.0','illuminaio':'0.52.0','IlluminaHumanMethylationEPICmanifest':'0.3.0'},
            'scientific_status':'CANNOT VERIFY','author_method_reproduction':False,'accepted_results':0,'model_calls':0,'new_clock_runs':0}
        self.json(self.out/'processing.json',self.processing)
        self.json(self.out/'worker-config.json',{'acquisition':str(self.acquisition_path),'output':str(self.out),'profile':self.profile})
        (self.out/'stdout.txt').write_text('Software fixture; no worker ran.\n');(self.out/'stderr.txt').write_text('')
        for role,name in {**self.outputs,'processing':'processing.json','worker_config':'worker-config.json','worker_stdout':'stdout.txt','worker_stderr':'stderr.txt'}.items():
            self.add(role,self.out/name)
        self.provenance={'schema':'clockbase-dnam-idat-provenance-v1','study_id':'GSE999001','sample_ids':self.ids,
            'profile':self.profile,'files':self.files,'elapsed_seconds':0,'scientific_status':'CANNOT VERIFY','accepted_results':0}
        self.provenance_path=self.out/'provenance.json';self.json(self.provenance_path,self.provenance)

    @staticmethod
    def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

    @staticmethod
    def json(path,value): Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

    @staticmethod
    def frame(path,frame,index=False):
        data=frame.to_csv(index=index,na_rep='NA',float_format='%.17g').encode()
        Path(path).write_bytes(gzip.compress(data,mtime=0) if str(path).endswith('.gz') else data)

    def add(self,role,path):
        path=Path(path).resolve();self.files[role]={'path':str(path),'sha256':self.sha(path),'size_bytes':path.stat().st_size}

    def reseal(self,role):
        self.add(role,self.files[role]['path']);self.json(self.provenance_path,self.provenance)

    def table(self,role): return pd.read_csv(self.files[role]['path'],index_col=0)

    def change_table(self,role,change):
        frame=self.table(role);change(frame);self.frame(self.files[role]['path'],frame,index=True);self.reseal(role)

    @staticmethod
    def binary(header,data):
        def string(text):
            value=text.encode('ascii');assert len(value)<128
            return bytes([len(value)])+value
        parts=[(1000,struct.pack('<I',len(data['address']))),(102,np.asarray(data['address'],dtype='<i4').tobytes()),
            (103,np.asarray(data['SD'],dtype='<u2').tobytes()),(104,np.asarray(data['Mean'],dtype='<u2').tobytes()),
            (107,np.asarray(data['NBeads'],dtype='u1').tobytes()),(400,struct.pack('<i',header['RedGreen'])),
            (402,string(header['Barcode'])),(403,string(header['ChipType'])),(404,string(header['position']))]
        offset=16+10*len(parts);directory=b'';body=b''
        for code,value in parts:
            directory+=struct.pack('<HQ',code,offset);body+=value;offset+=len(value)
        return b'IDAT'+struct.pack('<QI',3,len(parts))+directory+body

    def reseal_acquired_channel(self,suffix,binary,update_raw=True):
        compressed=Path(self.files['compressed:'+suffix]['path']);raw=Path(self.files['raw:'+suffix]['path'])
        compressed.write_bytes(gzip.compress(binary,mtime=0))
        if update_raw: raw.write_bytes(binary)
        receipt=Path(self.files['source_receipt:'+suffix]['path']);value=json.loads(receipt.read_text())
        value.update(sha256=self.sha(compressed),size_bytes=compressed.stat().st_size);self.json(receipt,value)
        self.json(self.files['decompression:'+suffix]['path'],{'compressed_sha256':self.sha(compressed),'sha256':self.sha(raw),'size_bytes':raw.stat().st_size})
        for kind in ('compressed','raw','source_receipt','decompression'):
            role=kind+':'+suffix;self.add(role,self.files[role]['path']);self.acquisition['files'][role]=copy.deepcopy(self.files[role])
        gsm,channel=suffix.split(':');pair=next(p for p in self.acquisition['pairs'] if p['sample_id']==gsm)
        pair['channels'][channel]['sha256']=self.sha(raw)
        self.acquisition['compressed_bytes']=sum(v['size_bytes'] for k,v in self.files.items() if k.startswith('compressed:'))
        self.json(self.acquisition_path,self.acquisition);self.reseal('acquisition')


class DnamIdatProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='disclosed-idat-software-fixture-')
        self.addCleanup(self.temp.cleanup);self.fixture=IdatFixture(self.temp.name)

    def verify(self): return verify_dnam_idat_provenance(self.fixture.provenance_path)

    def test_positive_independent_raw_mapping_arithmetic_and_mask_scope(self):
        value=self.verify();measurement=value['measurement']
        self.assertEqual(value['sample_ids'],IdatFixture.ids)
        self.assertEqual(set(value['files']),set(self.fixture.files)|{'provenance'})
        self.assertTrue(all(set(record)=={'path','sha256'} for record in value['files'].values()))
        self.assertEqual(measurement['beta_cells_compared'],6)
        self.assertEqual(measurement['masked_cells'],4)
        self.assertEqual(measurement['observed_cells'],2)
        self.assertEqual(measurement['normalization_validation'],'CANNOT VERIFY')
        self.assertEqual(measurement['specimen_identity_beyond_GEO'],'CANNOT VERIFY')
        self.assertEqual(measurement['accepted_results'],0)
        self.assertNotIn(str(self.fixture.root),json.dumps(measurement))
        self.assertEqual(verify_dnam_idat_provenance(self.fixture.provenance_path,self.fixture.files['matrix']['path']),value)

    def test_raw_redgreen_zero_is_retained_without_channel_inference(self):
        for channel in ('Grn','Red'):
            h,q=_idat_binary_readback(self.fixture.files['raw:'+IdatFixture.ids[0]+':'+channel]['path'])
            self.assertEqual(h['RedGreen'],0);self.assertEqual(len(q['address']),10)
        self.assertEqual(self.verify()['measurement']['status'],'PASS')

    def test_update_controller_accepts_exact_path_hash_contract(self):
        from src.update_pipeline_2026_09_20 import read_dnam_idat_provenance
        value=read_dnam_idat_provenance(self.fixture.provenance_path,
            self.fixture.sha(self.fixture.provenance_path),'GSE999001',
            sources={role:self.fixture.files[role] for role in ('matrix','series','samples')})
        self.assertEqual(value['scientific_status'],'CANNOT VERIFY')
        self.assertEqual(value['files']['wrapper'],value['files']['provenance'])

    def test_noncanonical_matrix_argument_rejected(self):
        other=self.fixture.root/'other.csv.gz';other.write_bytes(Path(self.fixture.files['matrix']['path']).read_bytes())
        with self.assertRaisesRegex(ValueError,'caller matrix'):verify_dnam_idat_provenance(self.fixture.provenance_path,other)

    def test_changed_implementation_hash_rejected(self):
        Path(self.fixture.files['implementation:rworker']['path']).write_text('changed')
        with self.assertRaisesRegex(ValueError,'source bytes changed'):self.verify()

    def test_missing_file_stops_with_path(self):
        path=Path(self.fixture.files['raw:'+IdatFixture.ids[0]+':Grn']['path']);path.unlink()
        with self.assertRaisesRegex(FileNotFoundError,'Cannot proceed - required data file not found'):self.verify()

    def test_omitted_channel_role_rejected(self):
        del self.fixture.files['raw:'+IdatFixture.ids[0]+':Red'];self.fixture.json(self.fixture.provenance_path,self.fixture.provenance)
        with self.assertRaisesRegex(ValueError,'acquisition'):self.verify()

    def test_unknown_provenance_role_rejected(self):
        p=self.fixture.root/'extra';p.write_text('software fixture');self.fixture.add('extra',p)
        self.fixture.json(self.fixture.provenance_path,self.fixture.provenance)
        with self.assertRaisesRegex(ValueError,'extra file roles'):self.verify()

    def test_source_url_cannot_move_between_gsm_owners(self):
        p=Path(self.fixture.files['samples']['path']);p.write_text(p.read_text().replace('/GSM100001/suppl/','/GSM100002/suppl/',1));self.fixture.reseal('samples')
        with self.assertRaisesRegex(ValueError,'exact GSM/channel ownership'):self.verify()

    def test_resealed_compressed_raw_mismatch_rejected(self):
        suffix=IdatFixture.ids[0]+':Grn';r=self.fixture.channel_data[suffix];data=copy.deepcopy(r['data']);data['Mean'][0]+=1
        self.fixture.reseal_acquired_channel(suffix,self.fixture.binary(r['header'],data),update_raw=False)
        with self.assertRaisesRegex(ValueError,'uncompressed bytes'):self.verify()

    def test_resealed_binary_intensity_disagrees_with_decoded_table(self):
        suffix=IdatFixture.ids[0]+':Grn';r=self.fixture.channel_data[suffix];data=copy.deepcopy(r['data']);data['Mean'][0]+=1
        self.fixture.reseal_acquired_channel(suffix,self.fixture.binary(r['header'],data))
        with self.assertRaisesRegex(ValueError,'decoded channel'):self.verify()

    def test_resealed_header_disagrees_with_binary(self):
        role='channel_header:'+IdatFixture.ids[0]+':Grn';p=Path(self.fixture.files[role]['path']);v=json.loads(p.read_text());v['position']='R08C01';self.fixture.json(p,v);self.fixture.reseal(role)
        with self.assertRaisesRegex(ValueError,'exported header'):self.verify()

    def test_swapped_raw_probe_values_rejected(self):
        def swap(frame):frame.iloc[[0,1],:]=frame.iloc[[1,0],:].to_numpy()
        self.fixture.change_table('raw_meth',swap)
        with self.assertRaisesRegex(ValueError,'binary manifest mappings'):self.verify()

    def test_swapped_sample_columns_rejected(self):
        self.fixture.change_table('raw_meth',lambda frame:setattr(frame,'columns',list(reversed(frame.columns))))
        with self.assertRaisesRegex(ValueError,'header/order'):self.verify()

    def test_swapped_probe_identities_rejected(self):
        self.fixture.change_table('meth_corrected',lambda frame:setattr(frame,'index',pd.Index(list(reversed(frame.index)),name=frame.index.name)))
        with self.assertRaisesRegex(ValueError,'probe/sample identities'):self.verify()

    def test_duplicate_physical_array_across_gsm_rejected(self):
        for channel in ('Grn','Red'):
            suffix=IdatFixture.ids[1]+':'+channel;r=self.fixture.channel_data[suffix]
            header={**r['header'],'position':'R01C01'}
            self.fixture.reseal_acquired_channel(suffix,self.fixture.binary(header,r['data']))
            role='channel_header:'+suffix
            self.fixture.json(self.fixture.files[role]['path'],{'sample_id':IdatFixture.ids[1],'channel':channel,**header})
            self.fixture.reseal(role)
        with self.assertRaisesRegex(ValueError,'duplicates physical arrays'):self.verify()

    def test_short_gsm_uses_geo_gsmnnn_directory(self):
        series=self.fixture.root/'ancient-series.soft';samples=self.fixture.root/'ancient-samples.soft'
        series.write_text('^SERIES = GSE1\n!Series_sample_id = GSM1\n')
        samples.write_text('^SAMPLE = GSM1\n!Sample_platform_id = GPL23976\n'+''.join(
            '!Sample_supplementary_file = https://ftp.ncbi.nlm.nih.gov/geo/samples/GSMnnn/GSM1/suppl/GSM1_Label_'+channel+'.idat.gz\n'
            for channel in ('Grn','Red')))
        self.assertEqual(_idat_literal_pairs('GSE1',series,samples)[0]['sample_id'],'GSM1')

    def test_resealed_beta_formula_change_rejected(self):
        self.fixture.change_table('beta_unmasked',lambda frame:frame.__setitem__(IdatFixture.ids[0],[0.6,0,np.nan]))
        with self.assertRaisesRegex(ValueError,'beta arithmetic'):self.verify()

    def test_exact_detection_threshold_is_not_masked(self):
        self.fixture.change_table('matrix',lambda frame:frame.__setitem__(IdatFixture.ids[0],[np.nan,np.nan,np.nan]))
        with self.assertRaisesRegex(ValueError,'masking differs'):self.verify()

    def test_failed_detection_cannot_be_unmasked(self):
        self.fixture.change_table('matrix',lambda frame:frame.__setitem__(IdatFixture.ids[1],[0.75,0.4,np.nan]))
        with self.assertRaisesRegex(ValueError,'masking differs'):self.verify()

    def test_probability_outside_unit_interval_rejected(self):
        self.fixture.change_table('detection_p',lambda frame:frame.__setitem__(IdatFixture.ids[0],[1.1,np.nan,0.001]))
        with self.assertRaisesRegex(ValueError,'probability range'):self.verify()

    def test_qc_count_mismatch_rejected(self):
        p=Path(self.fixture.files['sample_qc']['path']);f=pd.read_csv(p);f.loc[0,'masked_count']+=1;self.fixture.frame(p,f);self.fixture.reseal('sample_qc')
        with self.assertRaisesRegex(ValueError,'QC counts'):self.verify()

    def test_manifest_address_absent_from_binary_rejected(self):
        self.fixture.change_table('manifest_I',lambda frame:frame.__setitem__('AddressA',[999,3]))
        with self.assertRaisesRegex(ValueError,'address absent'):self.verify()

    def test_scientific_acceptance_cannot_be_promoted(self):
        self.fixture.provenance['scientific_status']='PASS';self.fixture.json(self.fixture.provenance_path,self.fixture.provenance)
        with self.assertRaisesRegex(ValueError,'scientific scope'):self.verify()

    def test_binary_directory_overlapping_offsets_rejected(self):
        p=Path(self.fixture.files['raw:'+IdatFixture.ids[0]+':Grn']['path']);data=bytearray(p.read_bytes());data[28:36]=data[18:26];p.write_bytes(data)
        with self.assertRaisesRegex(ValueError,'overlap'):_idat_binary_readback(p)


class SupportedIdatFixture(IdatFixture):
    """SOFTWARE FIXTURE: 7 full probes, 5 raw probes, 3 corrected probes.

    Two probes lack addresses and two further probes lack annotation groups.
    Values are deterministic integrity cases, not normalization measurements.
    """
    full_probes = ['cg00000001', 'cg00000002', 'cg00000004', 'cg00000005',
                   'cg00000003', 'cg00000006', 'cg00000007']
    raw_probes = ['cg00000001', 'cg00000002', 'cg00000004', 'cg00000003', 'cg00000006']

    def __init__(self, root):
        super().__init__(root)
        prior_hash = self.sha(self.acquisition_path)
        self.profile = _idat_verification_profile('minfi_ssnoob_beta100_detection01_supported_v2')
        self.plan['profile'] = self.profile
        self.json(self.plan_path, self.plan); self.add('protocol', self.plan_path)
        self.acquisition.update(profile=self.profile, plan_sha256=self.sha(self.plan_path),
                                reused_acquisition_sha256=prior_hash, new_downloads=0)
        self.json(self.acquisition_path, self.acquisition); self.add('acquisition', self.acquisition_path)
        self.processing['profile'] = self.profile; self.provenance['profile'] = self.profile
        self.json(self.files['worker_config']['path'],
                  {'acquisition': str(self.acquisition_path), 'output': str(self.out), 'profile': self.profile})
        self.add('worker_config', self.files['worker_config']['path'])
        self.write_table('manifest_I', pd.DataFrame([
            ['cg00000001', 1, 2, 'Grn'], ['cg00000002', 3, 4, 'Red'],
            ['cg00000004', 1, 2, 'Grn'], ['cg00000005', 11, 12, 'Red']],
            columns=['probe_id', 'AddressA', 'AddressB', 'Color']))
        self.write_table('manifest_II', pd.DataFrame([
            ['cg00000003', 5], ['cg00000006', 2], ['cg00000007', 13]], columns=['probe_id', 'AddressA']))
        self.new_output('annotation_groups', 'annotation-groups.csv', pd.DataFrame([
            ['cg00000001', 'I', 'Grn'], ['cg00000002', 'I', 'Red'],
            ['cg00000005', 'I', 'Red'], ['cg00000003', 'II', '']], columns=['probe_id', 'Type', 'Color']))
        self.new_output('probe_membership', 'probe-membership.csv', pd.DataFrame([
            ['cg00000001', 'I', 'Grn', 'TRUE', 'TRUE', 'retained'],
            ['cg00000002', 'I', 'Red', 'TRUE', 'TRUE', 'retained'],
            ['cg00000004', 'I', 'Grn', 'TRUE', 'FALSE', 'annotation_absent'],
            ['cg00000005', 'I', 'Red', 'FALSE', 'TRUE', 'address_absent'],
            ['cg00000003', 'II', '', 'TRUE', 'TRUE', 'retained'],
            ['cg00000006', 'II', '', 'TRUE', 'FALSE', 'annotation_absent'],
            ['cg00000007', 'II', '', 'FALSE', 'FALSE', 'address_absent']],
            columns=['probe_id', 'Type', 'Color', 'address_supported', 'annotation_supported', 'disposition']))
        self.new_output('missing_addresses', 'missing-addresses.csv', pd.DataFrame({'address': [11, 12, 13]}))
        for role, values in [('raw_meth', [[20, 120], [43, 143], [20, 120], [50, 150], [20, 120]]),
                             ('raw_unmeth', [[10, 110], [33, 133], [10, 110], [53, 153], [23, 123]])]:
            frame = pd.DataFrame(values, index=self.raw_probes, columns=self.ids); frame.index.name = 'probe_id'
            self.write_table(role, frame, index=True)
        raw_detection = pd.DataFrame([[0.01, 0.010001], [np.nan, 0.005], [0.02, 0.03],
                                      [0.001, np.inf], [np.nan, 0.0002]], index=self.raw_probes, columns=self.ids)
        raw_detection.index.name = 'probe_id'
        self.new_output('raw_detection_p', 'raw-detection-p.csv.gz', raw_detection, index=True)
        beta = self.table('beta_unmasked')
        self.new_output('beta_single_array', 'beta-single-array.csv.gz', beta, index=True)
        qc = pd.read_csv(self.files['sample_qc']['path'])
        qc['raw_median_M'] = [20, 120]; qc['raw_median_U'] = [23, 123]
        self.write_table('sample_qc', qc)
        self.processing['compatibility'].update(required_manifest_addresses=13, missing_manifest_addresses=3,
            typeI_count=4, typeII_count=3, raw_probe_count=5, annotation_excluded_probe_count=2,
            address_excluded_probe_count=2, corrected_probe_count=3)
        self.save_processing()

    def write_table(self, role, frame, index=False):
        self.frame(self.files[role]['path'], frame, index=index); self.reseal(role)

    def table(self, role):
        return pd.read_csv(self.files[role]['path'], index_col=0, float_precision='round_trip')

    def new_output(self, role, filename, frame, index=False):
        self.outputs[role] = filename; self.frame(self.out/filename, frame, index=index)
        self.add(role, self.out/filename)

    def save_processing(self):
        self.json(self.files['processing']['path'], self.processing); self.reseal('processing')


class SupportedDnamIdatProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='disclosed-idat-v2-software-fixture-')
        self.addCleanup(self.temp.cleanup); self.fixture = SupportedIdatFixture(self.temp.name)

    def verify(self): return verify_dnam_idat_provenance(self.fixture.provenance_path)

    def test_v2_separate_raw_and_corrected_rosters_and_path_free_contract(self):
        result = self.verify(); measured = result['measurement']
        self.assertTrue(all(set(value) == {'path', 'sha256'} for value in result['files'].values()))
        self.assertEqual(measured['full_manifest_probe_count'], 7)
        self.assertEqual(measured['raw_probe_count'], 5)
        self.assertEqual(measured['corrected_probe_count'], 3)
        self.assertEqual(measured['address_excluded_probe_count'], 2)
        self.assertEqual(measured['annotation_excluded_probe_count'], 2)
        self.assertEqual(measured['raw_methylation_cells_compared'], 20)
        self.assertEqual(measured['beta_cells_compared'], 6)
        self.assertEqual(measured['single_array_beta_cells_compared'], 6)
        self.assertEqual(measured['masked_cells'], 4)
        self.assertEqual(measured['sample_qc'][0]['raw_median_M'], 20)
        self.assertEqual(measured['sample_vs_joint_max_abs_delta'], dict.fromkeys(self.fixture.ids, 0.0))
        self.assertEqual(measured['acquisition_mode'], 'reused')
        self.assertEqual(measured['acquisition_reuse_scope'],
                         'Source bytes independently reread; prior acquisition receipt identified by hash only')
        for key in ['normalization_validation', 'manifest_annotation_validation', 'scientific_status']:
            self.assertEqual(measured[key], 'CANNOT VERIFY')
        self.assertNotIn(str(self.fixture.root), json.dumps(measured))
        self.assertEqual(measured['accepted_results'], 0)

    def test_v2_update_controller_accepts_path_hash_only_records(self):
        from src.update_pipeline_2026_09_20 import read_dnam_idat_provenance
        result = read_dnam_idat_provenance(self.fixture.provenance_path, self.fixture.sha(self.fixture.provenance_path),
            'GSE999001', sources={role: self.fixture.files[role] for role in ('matrix', 'series', 'samples')})
        self.assertEqual(result['measurement']['raw_probe_count'], 5)
        self.assertEqual(result['files']['wrapper'], result['files']['provenance'])

    def test_fresh_v2_acquisition_without_reuse_markers_is_supported(self):
        for name in ('reused_acquisition_sha256', 'new_downloads'):
            self.fixture.acquisition.pop(name)
        self.fixture.json(self.fixture.acquisition_path, self.fixture.acquisition); self.fixture.reseal('acquisition')
        measured = self.verify()['measurement']
        self.assertEqual(measured['acquisition_mode'], 'fresh')
        self.assertEqual(measured['acquisition_reuse_scope'],
                         'Fresh acquisition; source bytes independently reread, no prior receipt declared')
        self.assertEqual(measured['raw_probe_count'], 5)

    def test_partial_reuse_hash_without_download_count_is_rejected(self):
        self.fixture.acquisition.pop('new_downloads')
        self.fixture.json(self.fixture.acquisition_path, self.fixture.acquisition); self.fixture.reseal('acquisition')
        with self.assertRaisesRegex(ValueError, 'prior receipt hash'): self.verify()

    def test_partial_download_count_without_reuse_hash_is_rejected(self):
        self.fixture.acquisition.pop('reused_acquisition_sha256')
        self.fixture.json(self.fixture.acquisition_path, self.fixture.acquisition); self.fixture.reseal('acquisition')
        with self.assertRaisesRegex(ValueError, 'prior receipt hash'): self.verify()

    def test_v1_profile_is_not_silently_relaxed(self):
        fixture = self.fixture
        for value in [fixture.provenance, fixture.plan, fixture.acquisition, fixture.processing]:
            value['profile'] = _idat_verification_profile()
        fixture.json(fixture.plan_path, fixture.plan); fixture.add('protocol', fixture.plan_path)
        fixture.acquisition['plan_sha256'] = fixture.sha(fixture.plan_path)
        fixture.json(fixture.acquisition_path, fixture.acquisition); fixture.add('acquisition', fixture.acquisition_path)
        fixture.save_processing()
        with self.assertRaisesRegex(ValueError, 'extra file roles'): self.verify()

    def test_resealed_wrong_exclusion_reason_is_rejected(self):
        fixture = self.fixture
        frame = pd.read_csv(fixture.files['probe_membership']['path'], dtype=str, keep_default_na=False)
        frame.loc[2, 'disposition'] = 'address_absent'; fixture.write_table('probe_membership', frame)
        with self.assertRaisesRegex(ValueError, 'probe exclusions'): self.verify()

    def test_resealed_boolean_claim_cannot_create_address_support(self):
        fixture = self.fixture
        frame = pd.read_csv(fixture.files['probe_membership']['path'], dtype=str, keep_default_na=False)
        frame.loc[3, 'address_supported'] = 'TRUE'; fixture.write_table('probe_membership', frame)
        with self.assertRaisesRegex(ValueError, 'probe exclusions'): self.verify()

    def test_resealed_missing_address_list_must_be_complete_and_ordered(self):
        self.fixture.write_table('missing_addresses', pd.DataFrame({'address': [13, 11]}))
        with self.assertRaisesRegex(ValueError, 'missing address roster'): self.verify()

    def test_resealed_annotation_color_conflict_is_not_an_exclusion(self):
        frame = pd.read_csv(self.fixture.files['annotation_groups']['path'], keep_default_na=False)
        frame.loc[0, 'Color'] = 'Red'; self.fixture.write_table('annotation_groups', frame)
        with self.assertRaisesRegex(ValueError, 'groups conflict'): self.verify()

    def test_address_absent_probe_still_cannot_have_conflicting_annotation(self):
        frame = pd.read_csv(self.fixture.files['annotation_groups']['path'], keep_default_na=False)
        frame.loc[2, 'Color'] = 'Grn'; self.fixture.write_table('annotation_groups', frame)
        with self.assertRaisesRegex(ValueError, 'groups conflict'): self.verify()

    def test_duplicate_annotation_group_identity_rejected(self):
        frame = pd.read_csv(self.fixture.files['annotation_groups']['path'], keep_default_na=False)
        self.fixture.write_table('annotation_groups', pd.concat([frame, frame.iloc[[0]]], ignore_index=True))
        with self.assertRaisesRegex(ValueError, 'duplicate identities'): self.verify()

    def test_type_ii_annotation_color_must_be_empty(self):
        frame = pd.read_csv(self.fixture.files['annotation_groups']['path'], keep_default_na=False)
        frame.loc[3, 'Color'] = 'Grn'; self.fixture.write_table('annotation_groups', frame)
        with self.assertRaisesRegex(ValueError, 'invalid type/color'): self.verify()

    def test_missing_control_address_stops_supported_v2(self):
        frame = pd.read_csv(self.fixture.files['manifest_controls']['path'])
        frame.loc[0, 'Address'] = 999; self.fixture.write_table('manifest_controls', frame)
        with self.assertRaisesRegex(ValueError, 'every control address'): self.verify()

    def test_raw_matrix_keeps_annotation_excluded_values(self):
        frame = self.fixture.table('raw_meth').drop(index='cg00000004')
        self.fixture.write_table('raw_meth', frame, index=True)
        with self.assertRaisesRegex(ValueError, 'roster order'): self.verify()

    def test_resealed_raw_cells_on_annotation_excluded_probe_are_checked(self):
        self.fixture.change_table('raw_meth', lambda frame: frame.__setitem__(self.fixture.ids[0], [20, 43, 999, 50, 20]))
        with self.assertRaisesRegex(ValueError, 'binary manifest mappings'): self.verify()

    def test_swapped_raw_probe_rows_rejected_even_with_matching_values(self):
        frame = self.fixture.table('raw_meth').iloc[[0, 2, 1, 3, 4]]
        self.fixture.write_table('raw_meth', frame, index=True)
        with self.assertRaisesRegex(ValueError, 'roster order'): self.verify()

    def test_swapped_gsm_columns_rejected(self):
        self.fixture.change_table('raw_meth', lambda frame: setattr(frame, 'columns', list(reversed(frame.columns))))
        with self.assertRaisesRegex(ValueError, 'header/order'): self.verify()

    def test_annotation_missing_probe_cannot_be_falsely_corrected(self):
        for role in ['matrix', 'beta_unmasked', 'meth_corrected', 'unmeth_corrected', 'detection_p', 'beta_single_array']:
            frame = self.fixture.table(role)
            extra = frame.iloc[[0]].copy(); extra.index = pd.Index(['cg00000004'], name='probe_id')
            self.fixture.write_table(role, pd.concat([frame.iloc[:2], extra, frame.iloc[2:]]), index=True)
        with self.assertRaisesRegex(ValueError, 'roster order'): self.verify()

    def test_raw_to_selected_detection_mismatch_rejected(self):
        self.fixture.change_table('raw_detection_p', lambda frame: frame.__setitem__(self.fixture.ids[0], [0.001, np.nan, 0.02, 0.001, np.nan]))
        with self.assertRaisesRegex(ValueError, 'full-raw subset'): self.verify()

    def test_excluded_raw_detection_probability_is_still_checked(self):
        self.fixture.change_table('raw_detection_p', lambda frame: frame.__setitem__(self.fixture.ids[0], [0.01, np.nan, 1.5, 0.001, np.nan]))
        with self.assertRaisesRegex(ValueError, 'unit interval'): self.verify()

    def test_corrected_qc_cannot_use_raw_roster_denominator(self):
        frame = pd.read_csv(self.fixture.files['sample_qc']['path']); frame['probe_count'] = 5
        self.fixture.write_table('sample_qc', frame)
        with self.assertRaisesRegex(ValueError, 'QC counts'): self.verify()

    def test_raw_qc_median_cannot_use_corrected_roster(self):
        frame = pd.read_csv(self.fixture.files['sample_qc']['path']); frame['raw_median_M'] = [43, 143]
        self.fixture.write_table('sample_qc', frame)
        with self.assertRaisesRegex(ValueError, 'full-raw medians'): self.verify()

    def test_single_array_beta_values_are_compared(self):
        self.fixture.change_table('beta_single_array', lambda frame: frame.__setitem__(self.fixture.ids[0], [0.9, 0, np.nan]))
        with self.assertRaisesRegex(ValueError, 'beta delta'): self.verify()

    def test_single_array_nonfinite_pattern_cannot_change(self):
        self.fixture.change_table('beta_single_array', lambda frame: frame.__setitem__(self.fixture.ids[0], [0.5, 0, 0.5]))
        with self.assertRaisesRegex(ValueError, 'nonfinite identities'): self.verify()

    def test_recorded_delta_is_recomputed_instead_of_trusted(self):
        self.fixture.processing['sample_vs_joint_max_abs_delta'][self.fixture.ids[0]] = 5e-13
        self.fixture.save_processing()
        with self.assertRaisesRegex(ValueError, 'beta delta'): self.verify()

    def test_small_recorded_repeatability_delta_can_pass(self):
        frame = self.fixture.table('beta_single_array'); before = frame.iloc[0, 0]; frame.iloc[0, 0] += 5e-13
        delta = frame.iloc[0, 0]-before
        self.fixture.write_table('beta_single_array', frame, index=True)
        self.fixture.processing['sample_vs_joint_max_abs_delta'][self.fixture.ids[0]] = delta
        self.fixture.save_processing()
        self.assertAlmostEqual(self.verify()['measurement']['sample_vs_joint_max_abs_delta'][self.fixture.ids[0]], delta, places=16)

    def test_r_write_table_rounding_does_not_falsely_fail_saved_delta(self):
        frame = self.fixture.table('beta_single_array'); before = frame.iloc[1, 1]
        in_memory = before+9.995e-13
        frame.iloc[1, 1] = float(format(in_memory, '.15g'))
        self.fixture.write_table('beta_single_array', frame, index=True)
        self.fixture.processing['sample_vs_joint_max_abs_delta'][self.fixture.ids[1]] = in_memory-before
        self.fixture.save_processing()
        measured = self.verify()['measurement']
        self.assertLessEqual(abs(measured['sample_vs_joint_max_abs_delta'][self.fixture.ids[1]]-(in_memory-before)),
                             measured['single_array_delta_csv_rounding_bound'])

    def test_recorded_single_array_delta_above_frozen_limit_still_fails(self):
        self.fixture.processing['sample_vs_joint_max_abs_delta'][self.fixture.ids[0]] = 1.01e-12
        self.fixture.save_processing()
        with self.assertRaisesRegex(ValueError, 'Recorded single/joint'): self.verify()

    def test_reused_acquisition_digest_must_be_full_sha256(self):
        self.fixture.acquisition['reused_acquisition_sha256'] = 'not-a-hash'
        self.fixture.json(self.fixture.acquisition_path, self.fixture.acquisition); self.fixture.reseal('acquisition')
        with self.assertRaisesRegex(ValueError, 'prior receipt hash'): self.verify()

    def test_v2_reuse_cannot_claim_unbound_receipt_or_new_downloads(self):
        self.fixture.acquisition['new_downloads'] = 1
        self.fixture.json(self.fixture.acquisition_path, self.fixture.acquisition); self.fixture.reseal('acquisition')
        with self.assertRaisesRegex(ValueError, 'zero new downloads'): self.verify()

    def test_v2_compatibility_counts_close_both_exclusion_classes(self):
        self.fixture.processing['compatibility']['annotation_excluded_probe_count'] = 0
        self.fixture.save_processing()
        with self.assertRaisesRegex(ValueError, 'compatibility counts'): self.verify()

    def test_v2_compatibility_counts_are_integers(self):
        self.fixture.processing['compatibility']['raw_probe_count'] = 5.0
        self.fixture.save_processing()
        with self.assertRaisesRegex(ValueError, 'compatibility counts'): self.verify()


if __name__=='__main__': unittest.main()
