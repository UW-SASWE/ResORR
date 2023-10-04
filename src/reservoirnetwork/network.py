
import numpy as np
import pandas as pd
import xarray as xr

import networkx as nx

from .reservoir import Reservoir


class ReservoirNetwork(nx.DiGraph):
    def __init__(self, network, start_time, *args, **kwargs):
        super().__init__(network, *args, **kwargs)
        # self.data = xr.open_dataset(data) if isinstance(data, str) \
        #     else data if isinstance(data, xr.Dataset) \
        #     else ValueError("Must provide either path of data file or as xr.Dataset")
        self.data = xr.Dataset(
            coords={
                'node': list(self.nodes),
                'time': pd.date_range(start_time, periods=1, freq='1D')
            }
        ) 
        self.network = {node: Reservoir(node) for node in self.nodes}
        self.time = start_time

        self.FIRST_RUN = 1

    def create_field(self, var, fill_value=0.0):
        """Create a new field if not already present in self.data"""
        if var not in self.data.variables:
            self.data[var] = xr.DataArray(
                data=np.full((len(self.nodes), len(self.data.time)), fill_value),
                dims=['node', 'time'],
                coords={'node': self.nodes, 'time': self.data.time}
            )

    def insert_new_time_step(self, time):
        """Insert a new time step and fill with np.nan for all variables"""
        if time not in self.data.time:
            data_vars = {}
            for var in self.data.variables:
                if var not in ['node', 'time']:
                    data_vars[var] = (['node', 'time'], np.full((len(self.nodes), 1), np.nan))
            new_timestep_ds = xr.Dataset(data_vars=data_vars, coords={'node': self.nodes, 'time': [self.time]})
            # self.data = xr.merge([self.data, new_timestep_ds])
            self.data = xr.concat([self.data, new_timestep_ds], dim='time')
        else:
            raise ValueError(f"Time {time} already exists in data.")

    def update(self, forcings, dt=1, algorithm='simple', reservoir_algorithm='outlet'):
        """Update the reservoir network for one time step.

        Args:
            dt (int, optional): time step in days. Defaults to 1 day.
            algorithm (str, optional): Defaults to 'simple'.
                - hydraulic - outflow from reservoir is simulated to estimate storage change. 
                - simple - requires storage change and unregulated inflow of previous node.
        """

        if self.FIRST_RUN:
            self.create_field('inflow', np.nan)
            self.create_field('outflow', np.nan)
            self.create_field('regulated_runoff', np.nan)
            self.create_field('natural_runoff', np.nan)
            self.create_field('theoretical_natural_runoff', np.nan)
            self.create_field('storage', np.nan)
            self.create_field('storage_change', np.nan)
            self.create_field('regulation', np.nan)
            self.FIRST_RUN = 0
        else:
            self.time += pd.Timedelta(days=dt)
            self.insert_new_time_step(self.time)
        
        # insert forcings into data
        for var in forcings.variables:
            if var not in ['node', 'time']:
                # if var not in self.data.variables:
                self.data[var] = forcings[var]#.sel(time=self.time)
                # else:
                #     self.data[var].loc[dict(time=self.time)] = forcings[var].sel(time=self.time)
        
        if algorithm == 'hydraulic':
            self._alg_hydraulic(forcings, dt, reservoir_algorithm=reservoir_algorithm)

        if algorithm == 'hydraulic_travel_time':
            self._alg_hydraulic_travel_time(forcings, dt)

        if algorithm == 'simple':
            self._alg_simple(forcings, dt)
        
        if algorithm == 'simple_obs_outflow':
            self._alg_simple_obs_outflow(forcings, dt)

        if algorithm == 'simple_travel_time':
            self._alg_simple_travel_time(forcings, dt)

    def _alg_simple(self, forcings, dt):
        for node in list(nx.topological_sort(self)):
            storage_change = float(forcings['storage_change'].sel(node=node, time=self.time))
            theoretical_natural_runoff = float(forcings['theoretical_natural_runoff'].sel(node=node, time=self.time))

            self.data['theoretical_natural_runoff'].loc[dict(node=node, time=self.time)] = theoretical_natural_runoff

            upstream_dams = list(self.predecessors(node))
            natural_runoff = theoretical_natural_runoff
            regulated_runoff = 0.0
            if len(upstream_dams) > 0:
                regulated_runoff = sum([float(self.data['outflow'].sel(node=n, time=self.time)) for n in upstream_dams])
                natural_runoff -= sum([float(self.data['theoretical_natural_runoff'].sel(node=n, time=self.time)) for n in upstream_dams])

            inflow = max([0, float(natural_runoff + regulated_runoff)])
            outflow = max([0, inflow - storage_change])
            regulation = theoretical_natural_runoff - inflow

            self.data['inflow'].loc[dict(node=node, time=self.time)] = inflow
            self.data['outflow'].loc[dict(node=node, time=self.time)] = outflow
            self.data['regulation'].loc[dict(node=node, time=self.time)] = regulation
            self.data['natural_runoff'].loc[dict(node=node, time=self.time)] = natural_runoff
            self.data['regulated_runoff'].loc[dict(node=node, time=self.time)] = regulated_runoff
            self.data['storage_change'].loc[dict(node=node, time=self.time)] = storage_change
    
    def _alg_simple_obs_outflow(self, forcings, dt):
        for node in list(nx.topological_sort(self)):
            storage_change = float(forcings['storage_change'].sel(node=node, time=self.time))
            theoretical_natural_runoff = float(forcings['theoretical_natural_runoff'].sel(node=node, time=self.time))

            self.data['theoretical_natural_runoff'].loc[dict(node=node, time=self.time)] = theoretical_natural_runoff

            upstream_dams = list(self.predecessors(node))
            natural_runoff = theoretical_natural_runoff
            regulated_runoff = 0.0
            if len(upstream_dams) > 0:
                regulated_runoff = sum([float(self.data['obs_outflow'].sel(node=n, time=self.time)) for n in upstream_dams])
                natural_runoff -= sum([float(self.data['theoretical_natural_runoff'].sel(node=n, time=self.time)) for n in upstream_dams])

            inflow = max([0, float(natural_runoff + regulated_runoff)])
            outflow = max([0, inflow - storage_change])
            regulation = theoretical_natural_runoff - inflow

            self.data['inflow'].loc[dict(node=node, time=self.time)] = inflow
            self.data['outflow'].loc[dict(node=node, time=self.time)] = outflow
            self.data['regulation'].loc[dict(node=node, time=self.time)] = regulation
            self.data['natural_runoff'].loc[dict(node=node, time=self.time)] = natural_runoff
            self.data['regulated_runoff'].loc[dict(node=node, time=self.time)] = regulated_runoff
            self.data['storage_change'].loc[dict(node=node, time=self.time)] = storage_change

    def _alg_simple_travel_time(self, forcings, dt):
        for node in list(nx.topological_sort(self)):
            storage_change = float(forcings['storage_change'].sel(node=node, time=self.time))
            unregulated_inflow = float(forcings['unregulated_inflow'].sel(node=node, time=self.time))

            upstreams = list(self.predecessors(node))
            upstream_outflow = 0.0
            upstream_unregulated_inflow = 0.0
            if len(upstreams) > 0:
                time_lags = [self.time - pd.to_timedelta(round(self.get_edge_data(upstream, node)['travel_time']), 'd') for upstream in upstreams]
                upstream_outflow = sum([float(self.data['outflow'].sel(node=n, time=t)) for n, t in zip(upstreams, time_lags)])
                upstream_unregulated_inflow = sum([float(self.data['unregulated_inflow'].sel(node=n, time=t)) for n, t in zip(upstreams, time_lags)])
            
            regulated_inflow = max([0, float(unregulated_inflow - upstream_unregulated_inflow + upstream_outflow)])
            outflow = max([0, regulated_inflow - storage_change])

            self.data['regulated_inflow'].loc[dict(node=node, time=self.time)] = regulated_inflow
            self.data['outflow'].loc[dict(node=node, time=self.time)] = outflow

    def _alg_hydraulic(self, forcings, dt, reservoir_algorithm='outlet'):
        for node in list(nx.topological_sort(self)):
            # get reservoir object
            res = self.network[node]
            # get inflow for this node
            natural_runoff = (forcings['natural_runoff'].sel(node=node, time=slice(self.time - pd.Timedelta(days=dt), self.time)).mean() * dt).data # m3

            # if there are upstream nodes, the inflow 
            regulated_runoff = 0
            theoretical_natural_runoff = natural_runoff.copy()
            upstreams = list(self.predecessors(node))
            if len(upstreams) > 0:
                # sum the outflow from upstream nodes to the inflow of this node
                regulated_runoff += sum([self.data['outflow'].sel(node=n, time=self.time) for n in upstreams]).data
                theoretical_natural_runoff += sum([self.data['theoretical_natural_runoff'].sel(node=n, time=self.time) for n in upstreams]).data
            
            inflow = natural_runoff + regulated_runoff
            res.update(inflow, dt=1, algorithm=reservoir_algorithm) # run one time-step of the reservoir model

            res_data = res.dumpds()

            self.data['inflow'].loc[dict(node=node, time=self.time)] = res_data['inflow']
            self.data['outflow'].loc[dict(node=node, time=self.time)] = res_data['outflow']
            self.data['regulation'].loc[dict(node=node, time=self.time)] = theoretical_natural_runoff - res_data['inflow']
            self.data['storage_change'].loc[dict(node=node, time=self.time)] = res_data['storage_change']
            self.data['storage'].loc[dict(node=node, time=self.time)] = res_data['storage']
            self.data['regulated_runoff'].loc[dict(node=node, time=self.time)] = regulated_runoff
            self.data['theoretical_natural_runoff'].loc[dict(node=node, time=self.time)] = theoretical_natural_runoff

    def _alg_hydraulic_travel_time(self, forcings, dt):
        for node in list(nx.topological_sort(self)):
            # get reservoir object
            res = self.network[node]
            # get inflow for this node
            natural_runoff = (forcings['natural_runoff'].sel(node=node, time=slice(self.time - pd.Timedelta(days=dt), self.time)).mean() * dt).data # m3

            # if there are upstream nodes, the inflow 
            regulated_runoff = 0
            theoretical_natural_runoff = natural_runoff.copy()
            upstreams = list(self.predecessors(node))
            if len(upstreams) > 0:
                # sum the outflow from upstream nodes to the inflow of this node
                time_lags = [self.get_edge_data(upstream, node)['travel_time'] for upstream in upstreams]
                #### NOTE: method='nearest' is used as of now to get the nearest available upstream outflow value to the self.time-time_lag time.
                ####       this must be handled later, either by interpolating, or using np.nan for the first few time steps.
                regulated_runoff += sum([self.data['outflow'].sel(node=n, time=self.time-pd.Timedelta(lag, 'D')) for n, lag in zip(upstreams, time_lags)]).data
                theoretical_natural_runoff += sum([self.data['theoretical_natural_runoff'].sel(node=n, time=self.time-pd.Timedelta(lag, 'D')) for n, lag in zip(upstreams, time_lags)]).data
            
            inflow = natural_runoff + regulated_runoff
            res.update(inflow, 1) # run one time-step of the reservoir model
            res_data = res.dumpds()

            # update data
            self.data['inflow'].loc[dict(node=node, time=self.time)] = res_data['inflow']
            self.data['outflow'].loc[dict(node=node, time=self.time)] = res_data['outflow']
            self.data['regulation'].loc[dict(node=node, time=self.time)] = theoretical_natural_runoff - res_data['inflow']
            self.data['storage_change'].loc[dict(node=node, time=self.time)] = res_data['storage_change']
            self.data['storage'].loc[dict(node=node, time=self.time)] = res_data['storage']
            self.data['regulated_runoff'].loc[dict(node=node, time=self.time)] = regulated_runoff
            self.data['theoretical_natural_runoff'].loc[dict(node=node, time=self.time)] = theoretical_natural_runoff