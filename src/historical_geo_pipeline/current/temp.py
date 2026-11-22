import biolearn
import numpy as np
import pandas as pd
import requests
from biolearn.data_library import DataLibrary, GeoData
from biolearn.imputation import impute_from_standard
from biolearn.model_gallery import ModelGallery
from tqdm import tqdm


def calculate_clock_missing_percentage(methylation_site_report, clock_list, clock_models):
    # Extract the missing percentages from the report
    missing_percentages = methylation_site_report['missing']
    
    # Reshape missing percentages to a DataFrame with samples as columns
    missing_percentages_df = missing_percentages.to_frame().T
    
    # Dictionary to store missing percentage per clock
    clock_missing_percentages = {}
    
    for clock_name in clock_list:
        # Get the clock model
        clock_model = clock_models.get(clock_name)
        
        # Get the CpGs used by the clock
        clock_cpgs = clock_model.feature_names
        
        # Subset the missing percentage DataFrame to the clock's CpGs
        clock_missing_df = missing_percentages_df[clock_cpgs]
        
        # Calculate the mean missing percentage across the clock's CpGs
        mean_missing_percentage = clock_missing_df.mean(axis=1).values[0]
        
        # Store the result in the dictionary
        clock_missing_percentages[clock_name] = mean_missing_percentage
    
    # Convert the dictionary to a DataFrame
    clock_missing_percentages_df = pd.DataFrame.from_dict(clock_missing_percentages, orient='index', columns=['Missing Percentage'])
    
    return clock_missing_percentages_df
# Generate the methylation site report
model = ModelGallery().get("Horvathv1")
model.methylation_sites()

model = ModelGallery().get("GrimAgeV1")
model.methylation_sites()

report = data_quality_report(model.methylation_sites())
report.show()

# Calculate the missing percentage per clock
clock_missing_percentages_df = calculate_clock_missing_percentage(report, clock_list, clock_models)

# Print the results
print(clock_missing_percentages_df)
