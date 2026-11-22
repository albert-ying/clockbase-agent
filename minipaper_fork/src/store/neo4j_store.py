"""
Neo4j database interaction module
Implements CRUD operations for Neo4j database

Install dependencies:
pip install neo4j

Or add to pyproject.toml:
"neo4j>=5.0.0",
"""

from typing import Dict, List, Optional, Any, Union
from neo4j import GraphDatabase
from neo4j.exceptions import Neo4jError
from utils.logger import logger
from configs.global_config import NEO4J_CONFIG
from configs.config_cls import Neo4jConfig

class Neo4jStore:
    """Neo4j database interaction class providing full CRUD functionality"""
    
    def __init__(self, config: Neo4jConfig = None):
        """
        Initialize Neo4j connection
        
        Args:
            uri: Neo4j database URI, e.g. "bolt://localhost:7687"
            username: Username
            password: Password  
            database: Database name, defaults to "neo4j"
        """
        self.config = config or NEO4J_CONFIG
        self.uri = self.config.uri
        self.username = self.config.user
        self.password = self.config.password.get_secret_value()
        self.database = self.config.database
        self.driver = None
        
    def connect(self):
        """Establish database connection"""
        try:
            self.driver = GraphDatabase.driver(
                self.uri, 
                auth=(self.username, self.password)
            )
            logger.info(f"Successfully connected to Neo4j database: {self.uri}")
        except Neo4jError as e:
            logger.error(f"Failed to connect to Neo4j database: {e}")
            raise
            
    def close(self):
        """Close database connection"""
        if self.driver:
            self.driver.close()
            logger.info("Neo4j connection closed")
            
    def __enter__(self):
        """Context manager entry"""
        self.connect()
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit"""
        self.close()
        
    def execute_query(self, query: str, parameters: Dict = None) -> List[Dict]:
        """
        Execute a Cypher query
        
        Args:
            query: Cypher query statement
            parameters: Query parameters
            
        Returns:
            List of query results
        """
        if not self.driver:
            raise RuntimeError("Database not connected, please call connect() first")
            
        parameters = parameters or {}
        
        try:
            with self.driver.session(database=self.database) as session:
                result = session.run(query, parameters)
                return [record.data() for record in result]
        except Neo4jError as e:
            logger.error(f"Query execution failed: {query}, error: {e}")
            raise
            
    def execute_write_query(self, query: str, parameters: Dict = None) -> List[Dict]:
        """
        Execute a write query (CREATE, UPDATE, DELETE, etc.)
        
        Args:
            query: Cypher query statement
            parameters: Query parameters
            
        Returns:
            List of query results
        """
        if not self.driver:
            raise RuntimeError("Database not connected, please call connect() first")
            
        parameters = parameters or {}
        
        try:
            with self.driver.session(database=self.database) as session:
                result = session.write_transaction(lambda tx: tx.run(query, parameters))
                return [record.data() for record in result]
        except Neo4jError as e:
            logger.error(f"Write query execution failed: {query}, error: {e}")
            raise

    # ============ CREATE Operations ============
    
    def create_node(self, label: str, properties: Dict[str, Any]) -> Dict:
        """
        Create a node
        
        Args:
            label: Node label
            properties: Node properties dictionary
            
        Returns:
            Created node information
        """
        # Build property string
        prop_str = ", ".join([f"{key}: ${key}" for key in properties.keys()])
        query = f"CREATE (n:{label} {{{prop_str}}}) RETURN n"
        
        result = self.execute_write_query(query, properties)
        if result:
            return result[0]['n']
        return {}
        
    def create_relationship(self, 
                          from_node_label: str, 
                          from_node_props: Dict[str, Any],
                          to_node_label: str, 
                          to_node_props: Dict[str, Any],
                          relationship_type: str,
                          relationship_props: Dict[str, Any] = None) -> Dict:
        """
        Create a relationship
        
        Args:
            from_node_label: Source node label
            from_node_props: Source node match properties
            to_node_label: Target node label  
            to_node_props: Target node match properties
            relationship_type: Relationship type
            relationship_props: Relationship properties
            
        Returns:
            Created relationship information
        """
        relationship_props = relationship_props or {}
        
        # Build match conditions
        from_match = ", ".join([f"{key}: ${key}" for key in from_node_props.keys()])
        to_match = ", ".join([f"{key}: $to_{key}" for key in to_node_props.keys()])
        
        # Build relationship properties
        rel_props_str = ""
        if relationship_props:
            rel_prop_str = ", ".join([f"{key}: $rel_{key}" for key in relationship_props.keys()])
            rel_props_str = f" {{{rel_prop_str}}}"
        
        query = f"""
        MATCH (from:{from_node_label} {{{from_match}}}), 
              (to:{to_node_label} {{{to_match}}})
        CREATE (from)-[r:{relationship_type}{rel_props_str}]->(to)
        RETURN r
        """
        
        # Prepare parameters
        parameters = from_node_props.copy()
        parameters.update({f"to_{k}": v for k, v in to_node_props.items()})
        parameters.update({f"rel_{k}": v for k, v in relationship_props.items()})
        
        result = self.execute_write_query(query, parameters)
        if result:
            return result[0]['r']
        return {}

    # ============ READ Operations ============
    
    def find_nodes(self, label: str = None, properties: Dict[str, Any] = None, 
                   limit: int = 100) -> List[Dict]:
        """
        Find nodes
        
        Args:
            label: Node label, optional
            properties: Match properties, optional
            limit: Result limit count
            
        Returns:
            List of matching nodes
        """
        properties = properties or {}
        
        # Build query
        if label:
            match_clause = f"MATCH (n:{label}"
        else:
            match_clause = "MATCH (n"
            
        if properties:
            prop_str = ", ".join([f"{key}: ${key}" for key in properties.keys()])
            match_clause += f" {{{prop_str}}}"
            
        match_clause += ")"
        
        query = f"{match_clause} RETURN n LIMIT {limit}"
        
        result = self.execute_query(query, properties)
        return [record['n'] for record in result]
        
    def find_node_by_id(self, node_id: int) -> Optional[Dict]:
        """
        Find a node by ID
        
        Args:
            node_id: Node ID
            
        Returns:
            Node information or None
        """
        query = "MATCH (n) WHERE id(n) = $node_id RETURN n"
        result = self.execute_query(query, {"node_id": node_id})
        
        if result:
            return result[0]['n']
        return None
        
    def find_relationships(self, relationship_type: str = None, 
                         from_label: str = None, to_label: str = None,
                         limit: int = 100) -> List[Dict]:
        """
        Find relationships
        
        Args:
            relationship_type: Relationship type, optional
            from_label: Source node label, optional
            to_label: Target node label, optional
            limit: Result limit count
            
        Returns:
            List of matching relationships
        """
        # Build match pattern
        from_pattern = f"(from:{from_label})" if from_label else "(from)"
        to_pattern = f"(to:{to_label})" if to_label else "(to)"
        
        if relationship_type:
            rel_pattern = f"-[r:{relationship_type}]->"
        else:
            rel_pattern = "-[r]->"
            
        query = f"""
        MATCH {from_pattern}{rel_pattern}{to_pattern}
        RETURN from, r, to
        LIMIT {limit}
        """
        
        return self.execute_query(query)
        
    def get_node_relationships(self, node_label: str, node_props: Dict[str, Any],
                             direction: str = "both") -> List[Dict]:
        """
        Get all relationships of a node
        
        Args:
            node_label: Node label
            node_props: Node match properties
            direction: Relationship direction ("incoming", "outgoing", "both")
            
        Returns:
            List of relationships
        """
        prop_str = ", ".join([f"{key}: ${key}" for key in node_props.keys()])
        
        if direction == "incoming":
            pattern = f"(other)-[r]->(n:{node_label} {{{prop_str}}})"
        elif direction == "outgoing":  
            pattern = f"(n:{node_label} {{{prop_str}}})-[r]->(other)"
        else:  # both
            pattern = f"(n:{node_label} {{{prop_str}}})-[r]-(other)"
            
        query = f"MATCH {pattern} RETURN r, other"
        
        return self.execute_query(query, node_props)

    # ============ UPDATE Operations ============
    
    def update_node(self, label: str, match_props: Dict[str, Any], 
                   update_props: Dict[str, Any]) -> List[Dict]:
        """
        Update node properties
        
        Args:
            label: Node label
            match_props: Match condition properties
            update_props: Properties to update
            
        Returns:
            List of updated nodes
        """
        # Build match conditions
        match_str = ", ".join([f"{key}: $match_{key}" for key in match_props.keys()])
        
        # Build update statement
        update_str = ", ".join([f"n.{key} = $update_{key}" for key in update_props.keys()])
        
        query = f"""
        MATCH (n:{label} {{{match_str}}})
        SET {update_str}
        RETURN n
        """
        
        # Prepare parameters
        parameters = {f"match_{k}": v for k, v in match_props.items()}
        parameters.update({f"update_{k}": v for k, v in update_props.items()})
        
        result = self.execute_write_query(query, parameters)
        return [record['n'] for record in result]
        
    def update_relationship(self, from_label: str, from_props: Dict[str, Any],
                          to_label: str, to_props: Dict[str, Any],
                          rel_type: str, update_props: Dict[str, Any]) -> List[Dict]:
        """
        Update relationship properties
        
        Args:
            from_label: Source node label
            from_props: Source node match properties
            to_label: Target node label
            to_props: Target node match properties
            rel_type: Relationship type
            update_props: Relationship properties to update
            
        Returns:
            List of updated relationships
        """
        # Build match conditions
        from_match = ", ".join([f"{key}: $from_{key}" for key in from_props.keys()])
        to_match = ", ".join([f"{key}: $to_{key}" for key in to_props.keys()])
        
        # Build update statement
        update_str = ", ".join([f"r.{key} = $update_{key}" for key in update_props.keys()])
        
        query = f"""
        MATCH (from:{from_label} {{{from_match}}})-[r:{rel_type}]->(to:{to_label} {{{to_match}}})
        SET {update_str}
        RETURN r
        """
        
        # Prepare parameters
        parameters = {f"from_{k}": v for k, v in from_props.items()}
        parameters.update({f"to_{k}": v for k, v in to_props.items()})
        parameters.update({f"update_{k}": v for k, v in update_props.items()})
        
        result = self.execute_write_query(query, parameters)
        return [record['r'] for record in result]

    # ============ DELETE Operations ============
    
    def delete_node(self, label: str, properties: Dict[str, Any], 
                   force: bool = False) -> int:
        """
        Delete a node
        
        Args:
            label: Node label
            properties: Match properties
            force: Whether to force delete (also delete relationships)
            
        Returns:
            Number of deleted nodes
        """
        prop_str = ", ".join([f"{key}: ${key}" for key in properties.keys()])
        
        if force:
            # Delete node along with its relationships
            query = f"MATCH (n:{label} {{{prop_str}}}) DETACH DELETE n"
        else:
            # Delete node only (node must have no relationships)
            query = f"MATCH (n:{label} {{{prop_str}}}) DELETE n"
        
        try:
            self.execute_write_query(query, properties)
            # Get deleted count
            count_query = f"MATCH (n:{label} {{{prop_str}}}) RETURN count(n) as count"
            result = self.execute_query(count_query, properties)
            return result[0]['count'] if result else 0
        except Neo4jError as e:
            if "Cannot delete node" in str(e) and not force:
                raise ValueError("Node has existing relationships and cannot be deleted. Set force=True to force delete or remove relationships first.")
            raise
            
    def delete_relationship(self, from_label: str, from_props: Dict[str, Any],
                          to_label: str, to_props: Dict[str, Any],
                          rel_type: str) -> int:
        """
        Delete a relationship
        
        Args:
            from_label: Source node label
            from_props: Source node match properties
            to_label: Target node label
            to_props: Target node match properties
            rel_type: Relationship type
            
        Returns:
            Number of deleted relationships
        """
        from_match = ", ".join([f"{key}: $from_{key}" for key in from_props.keys()])
        to_match = ", ".join([f"{key}: $to_{key}" for key in to_props.keys()])
        
        query = f"""
        MATCH (from:{from_label} {{{from_match}}})-[r:{rel_type}]->(to:{to_label} {{{to_match}}})
        DELETE r
        RETURN count(r) as count
        """
        
        # Prepare parameters
        parameters = {f"from_{k}": v for k, v in from_props.items()}
        parameters.update({f"to_{k}": v for k, v in to_props.items()})
        
        result = self.execute_write_query(query, parameters)
        return result[0]['count'] if result else 0
        
    def clear_database(self, confirm: bool = False):
        """
        Clear the entire database (dangerous operation)
        
        Args:
            confirm: Confirm the operation
        """
        if not confirm:
            raise ValueError("Clearing the database requires setting confirm=True to confirm the operation")
            
        query = "MATCH (n) DETACH DELETE n"
        self.execute_write_query(query)
        logger.warning("Database has been cleared")

    # ============ Utility Methods ============
    
    def get_database_info(self) -> Dict[str, Any]:
        """Get basic database information"""
        queries = {
            "node_count": "MATCH (n) RETURN count(n) as count",
            "relationship_count": "MATCH ()-[r]->() RETURN count(r) as count",
            "labels": "CALL db.labels() YIELD label RETURN collect(label) as labels",
            "relationship_types": "CALL db.relationshipTypes() YIELD relationshipType RETURN collect(relationshipType) as types"
        }
        
        info = {}
        for key, query in queries.items():
            try:
                result = self.execute_query(query)
                if result:
                    if key in ["labels", "relationship_types"]:
                        info[key] = result[0][key.split('_')[0] if key == "relationship_types" else key]
                    else:
                        info[key] = result[0]["count"]
            except Exception as e:
                info[key] = f"Failed to retrieve: {e}"
                
        return info
        
    def test_connection(self) -> bool:
        """Test database connection"""
        try:
            result = self.execute_query("RETURN 1 as test")
            return result[0]["test"] == 1
        except Exception:
            return False


# ============ Usage Example ============

if __name__ == "__main__":
       
    # Use context manager
    with Neo4jStore() as store:
        # Test connection
        if store.test_connection():
            print("✅ Neo4j connection successful")
        else:
            print("❌ Neo4j connection failed")
            exit(1)
            
        
        # Query nodes
        paper = store.find_nodes("Paper", {"title": "Lactobacillus plantarum TWK10 Attenuates Aging-Associated Muscle Weakness, Bone Loss, and Cognitive Impairment by Modulating the Gut Microbiome in Mice"})
        print(f"Paper: {paper}")
        pdf_url = paper[0]['pdf_url']
        doi = paper[0]['doi']
        print(f"pdf_url: {pdf_url}")
        print(f"doi: {doi}")
        

        
        
