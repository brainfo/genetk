#!/usr/bin/env python
"""
Test script for GONetwork class functionality.
"""

import os
import sys
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for testing

def test_go_network():
    """Test GONetwork functionality with the example data."""
    
    # Import the GONetwork class
    try:
        from ..src.genetk.network.go_network import GONetwork
        print("✓ Successfully imported GONetwork")
    except ImportError as e:
        print(f"✗ Failed to import GONetwork: {e}")
        return False
    
    # Test data file path
    test_file = "GO_Biological_Process_2023.human.enrichr.reports.txt"
    
    if not os.path.exists(test_file):
        print(f"✗ Test data file not found: {test_file}")
        return False
    
    print(f"✓ Found test data file: {test_file}")
    
    try:
        # Initialize network
        print("Testing GONetwork initialization...")
        go_net = GONetwork(test_file)
        print(f"✓ Network initialized with {go_net.G.number_of_nodes()} nodes and {go_net.G.number_of_edges()} edges")
        
        # Test basic network properties
        print(f"✓ Data shape: {go_net.data.shape}")
        print(f"✓ Sample node attributes: {list(go_net.G.nodes(data=True))[0]}")
        print(f"✓ Sample edge attributes: {list(go_net.G.edges(data=True))[0] if go_net.G.edges() else 'No edges'}")
        
        # Test filtering by node attributes
        print("\nTesting node filtering...")
        filtered_net = go_net.filter_nodes_by_attribute('odds_ratio')
        print(f"✓ Filtered network: {filtered_net.G.number_of_nodes()} nodes (from {go_net.G.number_of_nodes()})")
        
        # Test filtering by edge weights
        print("\nTesting edge filtering...")
        edge_filtered_net = filtered_net.filter_edges_by_weight()
        print(f"✓ Edge filtered network: {edge_filtered_net.G.number_of_edges()} edges (from {filtered_net.G.number_of_edges()})")
        
        # Test community detection
        print("\nTesting community detection...")
        communities = edge_filtered_net.community()
        num_communities = len(set(communities.values()))
        print(f"✓ Found {num_communities} communities")
        
        # Create figures directory if it doesn't exist
        os.makedirs('figures', exist_ok=True)
        
        # Test visualization with community coloring
        print("\nTesting visualization with community coloring...")
        edge_filtered_net.bubble(
            name='test_community',
            community_colored=True,
            text_annot=True
        )
        print("✓ Community-colored visualization created")
        
        # Test visualization with p-value coloring
        print("\nTesting visualization with p-value coloring...")
        edge_filtered_net.bubble(
            name='test_pvalue',
            text_annot=True,
            color_attribute='adjusted_p_value',
            color='red'
        )
        print("✓ P-value colored visualization created")
        
        # Test with specific term list
        print("\nTesting with specific term list...")
        # Get a few terms from the network
        sample_terms = list(edge_filtered_net.G.nodes())[:3]
        edge_filtered_net.bubble(
            name='test_terms',
            text_annot=True,
            term_list=sample_terms,
            color='blue'
        )
        print(f"✓ Visualization with specific terms created: {sample_terms}")
        
        print("\n🎉 All tests passed successfully!")
        return True
        
    except Exception as e:
        print(f"✗ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_go_network()
    sys.exit(0 if success else 1)