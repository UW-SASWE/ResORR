import numpy as np
import numpy.ma as ma
import pandas as pd

from landlab import Component, FieldError
from landlab.components import FlowDirectorSteepest
from landlab.grid.mappers import map_link_head_node_to_link

class StreamflowRegulation(Component):
    """_brief_summary_
    _detailed_summary_
    Construction:
        
        Component (_type_): _description_
    Parameters
    ----------
        fluxes: DataRecord that has modeled total inflow and storage change for each reservoir in 
                the network
    """
    _name = "StreamflowRegulation"

    _unit_agnostic = False

    _info = {
        "reservoir__total_inflow": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m3/d",
            "mapping": "node"
        },
        "reservoir__storage_change": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m3/d",
            "mapping": "node"
        },
        "reservoir__release": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m3/d",
            "mapping": "node"
        },
        "reservoir__regulated_inflow": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m3/d",
            "mapping": "node"
        },
        "reservoir__unregulated_inflow": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m3/d",
            "mapping": "node"
        },
        "reservoir__abstract_elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node"
        },
        "river__storage": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m3",
            "mapping": "link"
        },
        "river__hydraulic_radius": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "link"
        },
        "river__release_inflow": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m3/d",
            "mapping": "link"
        },
        "river__slope": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/m",
            "mapping": "link"
        },
        "river__length": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "link"
        },
        "river__roughness": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "link"
        }
    }

    def __init__(
        self, grid#, fluxes
    ):
        """_summary_
        Args:
            grid (_type_): _description_
            start_time (_type_): _description_
        """
        super().__init__(grid)

        # self._grid = grid

        self.initialize_output_fields()
        
        # self._fluxes = fluxes

        self._time = 0
        self._time_idx = 0

        self._inflow = grid.at_node["reservoir__total_inflow"]
        self._abstract_elevation = grid.at_node["reservoir__abstract_elevation"]
        self._storage_change = grid.at_node["reservoir__storage_change"]

        self._time_variable_flux_attributes = [
            "modeled_inflow",
            "storage_change"
        ]


    def _calc_outflow(self):
        """Calculates outflow using mass balance: O = I - ∆S
        """
        self._outflow = self._grid.at_node['reservoir__total_inflow'] - self._grid.at_node['reservoir__storage_change']
        self._grid.at_node['reservoir__release'] = self._outflow
        
        return self._outflow

    # def _create_new_time(self):
    #     """_summary_
    #     """
    #     if self._time_idx != 0:
    #         self._fluxes.add_record(time=[self._time])
            
    #         self._fluxes.ffill_grid_element_and_id()

    #         # copy parcel attributes forward in time.
    #         for at in self._time_variable_flux_attributes:
    #             self._fluxes.dataset[at].values[
    #                 :, self._time_idx
    #             ] = self._fluxes.dataset[at].values[:, self._time_idx - 1]

    #     self._this_timesteps_fluxes = np.zeros_like(
    #         self._fluxes.dataset.element_id, dtype=bool
    #     )
    #     self._this_timesteps_fluxes[:, -1] = True


    def _calc_regulated_inflow(self):
        """Calculates the upstream inflow
        """
        # flow director needed for finding upstream nodes. Use the abstract_elevation field which
        # can have any abstract value, with the end node having the lowest abstract_elevation, and
        # increasing upstream
        fdr = FlowDirectorSteepest(self._grid, self._abstract_elevation)
        fdr.run_one_step()

        # find the links that contribute to each node
        upstream_contributing_links_at_node = np.where(
            fdr.flow_link_incoming_at_node() == 1, self._grid.links_at_node, -1
        )

        # # find upstream links
        # upstream_links = ma.MaskedArray(
        #     upstream_contributing_links_at_node, 
        #     mask=upstream_contributing_links_at_node==-1
        # )
        # # find upstream nodes
        # upstream_nodes = ma.MaskedArray(
        #     fdr.upstream_node_at_link()[upstream_links], 
        #     mask=upstream_links.mask
        # )
        
        # calculate regulated inflow as the sum of release for now, later use links and link 
        # storage to calculate regulated inflow
        self._grid.at_link['river__regulated_flow']
        regulated_inflow = ma.sum(
            ma.MaskedArray(
                self._grid.at_link['river__regulated_flow'].reshape(-1, 1)[upstream_contributing_links_at_node],
                mask = upstream_contributing_links_at_node == -1
            ), axis=1 
        ).filled(0)

        self._grid.at_node['reservoir__regulated_inflow'] = regulated_inflow

        return regulated_inflow

    def _calc_unregulated_inflow(self):
        """Calculates the unregulated inflow reaching a reservoir
        """
        total_inflow = self._grid.at_node['reservoir__total_inflow']
        regulated_inflow = self._grid.at_node['reservoir__regulated_inflow']

        unregulated_inflow = total_inflow - regulated_inflow
        self._grid.at_node['reservoir__unregulated_inflow'] = unregulated_inflow

        return unregulated_inflow

    def _calc_k(self):
        return np.power(self._grid.at_link["river__hydraulic_radius"], 2/3) * np.sqrt(self._grid.at_link["river__slope"])/(self._grid.at_link["river__roughness"]*self._grid.at_link["river__length"])
    
    def _route_through_stream(self):
        """Calculate inflow to river
        """
        # Qin = self._grid.at_node["reservoir__release"][self._grid.node_at_link_head])
        Qin = map_link_head_node_to_link(self._grid, "reservoir__release")
        self._grid.at_link["river__release_inflow"] = Qin
        
        V_t = self._grid.at_link["river__storage"]

        k = self._calc_k()

        V_t_1 = (Qin/k) + (V_t - Qin/k) * np.exp(-k) # ∆t = 1
        self._grid.at_link["river__storage"] = V_t_1

        # calculate outflow to next node
        Qout = Qin - (V_t_1 - V_t)
        self.grid.at_link["river__regulated_flow"] = Qout

        return Qout
    
    def _calc_upstream_storage_change(self):
        """Calculate cumulative upstream storage change
        """
        # flow director needed for finding upstream nodes. Use the abstract_elevation field which
        # can have any abstract value, with the end node having the lowest abstract_elevation, and
        # increasing upstream
        fdr = FlowDirectorSteepest(self._grid, self._abstract_elevation)
        fdr.run_one_step()

        # find the links that contribute to each node
        upstream_contributing_links_at_node = np.where(
            fdr.flow_link_incoming_at_node() == 1, self._grid.links_at_node, -1
        )

        # find upstream links
        upstream_links = ma.MaskedArray(
            upstream_contributing_links_at_node, 
            mask=upstream_contributing_links_at_node==-1
        )
        # find upstream nodes
        upstream_nodes = ma.MaskedArray(
            fdr.upstream_node_at_link()[upstream_links], 
            mask=upstream_links.mask
        )
        
        # calculate regulated inflow as the sum of release for now, later use links and link 
        # storage to calculate regulated inflow
        upstream_storage_change = ma.sum(
            ma.MaskedArray(
                self._grid.at_node['reservoir__storage_change'].reshape(-1, 1)[upstream_nodes],
                mask = upstream_nodes.mask
            ), axis=1 
        ).filled(0)

        self._grid.at_node['reservoir__upstream_storage_change'] = upstream_storage_change

        return upstream_storage_change


    def run_one_step(self, time, inflow, storage_change):
        """_summary_
        Args:
            dt (int, optional): _description_. Defaults to 1.
        """
        self._grid.at_node["reservoir__total_inflow"] = inflow
        self._grid.at_node["reservoir__storage_change"] = storage_change

        if time is not None:
            self._time = time
        else:
            self._time = self._time + pd.Timedelta(1, 'D')
        self._time_idx += 1

        self._calc_outflow()
        self._route_through_stream()
        self._calc_regulated_inflow()
        self._calc_unregulated_inflow()
        self._calc_upstream_storage_change()