import pytest
from pathlib import Path
import xarray as xr
import pandas as pd


def test_import():
    import resorr

# test cumberland
@pytest.fixture(scope="session")
def cumberland_save_dir(tmp_path_factory):
    return tmp_path_factory.mktemp('cumberland')

@pytest.fixture(scope="session")
def generate_network(cumberland_save_dir):
    ### NETWORK GENERATION
    flow_dir_fn = 'tests/data/cumberland/flow-directions.tif'
    stations_fn = 'tests/data/cumberland/cumberland-dams.csv'

    from resorr.generate_network import generate_network

    G = generate_network(flow_dir_fn, stations_fn, cumberland_save_dir)

    return G

@pytest.fixture(scope="session")
def prepare_input_data(generate_network, cumberland_save_dir):
    ### DATA PREP
    from resorr.data_prep import generate_forcings_from_rat

    inflow_dir = 'tests/data/cumberland/unregulated-inflow'
    storage_change_dir = 'tests/data/cumberland/satellite-dels'

    rat_output_level = 'rat_outputs'

    resorr_forcings = generate_forcings_from_rat(
        generate_network, inflow_dir, storage_change_dir, 
        cumberland_save_dir, rat_output_level=rat_output_level
    )

    return resorr_forcings

@pytest.fixture(scope="session")
def run_resorr(cumberland_save_dir, generate_network, prepare_input_data):
    from resorr.network import ReservoirNetwork

    start_time = pd.to_datetime('2019-01-01')
    end_time = pd.to_datetime('2019-01-31')

    forcings = prepare_input_data.sel(time=slice(start_time, end_time))
    G = generate_network

    reservoir_network = ReservoirNetwork(G, start_time)

    for timestep in forcings.time.values:
        dt = forcings['dt'].sel(time=timestep).values.item()
        reservoir_network.update(forcings, dt, 'wb')

    return reservoir_network


def test_cumberland(cumberland_save_dir, generate_network):
    G = generate_network

    edges = G.edges_to_gdf()
    nodes = G.nodes_to_gdf()
    
    # for cumberland, there should be 7 edges, 8 nodes
    assert len(edges) == 7
    assert len(nodes) == 8
    # check if files were created
    pts_save_fn = cumberland_save_dir / 'rivreg_network_pts.shp'
    save_fn = cumberland_save_dir / 'rivreg_network.shp'
    assert Path(pts_save_fn).exists()
    assert Path(save_fn).exists()

def test_cumberland_data_prep(cumberland_save_dir, prepare_input_data):
    # check if files were created
    regulation_inp_data_fn = cumberland_save_dir / 'resorr_forcings.nc'
    assert Path(regulation_inp_data_fn).exists()

    resorr_forcings = prepare_input_data
    # check if TNR and stprage change are present
    assert 'theoretical_natural_runoff' in resorr_forcings
    assert 'storage_change' in resorr_forcings

def test_resorr_run(run_resorr):
    resorr_output = run_resorr.data

    # check if all variables are present
    assert 'inflow' in resorr_output
    assert 'outflow' in resorr_output
    assert 'regulated_runoff' in resorr_output
    assert 'natural_runoff' in resorr_output
    assert 'regulation' in resorr_output