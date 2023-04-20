import re

import h5py
import numpy as np
import fsspec
import pandas as pd
import anndata

from cirrocumulus.ot.population import Population
from cirrocumulus.ot.util import chain_transport_maps, find_path, unique_timepoint
from cirrocumulus.util import get_scheme


class TransportMapModel:
    def __init__(self, tmaps, meta, timepoints=None, day_pairs=None, cache=False):
        self.tmaps = tmaps
        self.meta = meta
        self.cache = cache
        if timepoints is None:
            timepoints = sorted(meta["day"].unique())
        self.timepoints = timepoints

        if day_pairs is None:
            day_pairs = [(timepoints[i], timepoints[i + 1]) for i in range(len(timepoints) - 1)]
        self.day_pairs = day_pairs

    def trajectories(self, populations):
        """Computes a trajectory for each population.

        Parameters
        ----------
        self : wot.TransportMapModel
            The TransportMapModel used to find ancestors and descendants of the population
        populations : list of wot.Population
            The target populations such as ones from self.population_from_cell_sets. THe populations must be from the same time.

        Returns
        -------
        trajectories : anndata.AnnData
            Rows : all cells, Columns : populations index. At point (i, j) : the probability that cell i is an
            ancestor/descendant of population j
        """
        unique_timepoint(*populations)  # check for unique timepoint
        trajectories = []

        populations = Population.copy(*populations, normalize=True, add_missing=False)
        population_names = [p.name for p in populations]
        initial_populations = populations

        def update(head, populations_to_update):
            idx = 0 if head else len(trajectories)
            trajectories.insert(idx, np.array([pop.p for pop in populations_to_update]).T)

        update(True, populations)
        while self.can_pull_back(*populations):
            populations = self.pull_back(*populations, as_list=True)
            update(True, populations)
        populations = initial_populations
        while self.can_push_forward(*populations):
            populations = self.push_forward(*populations, as_list=True)
            update(False, populations)

        return anndata.AnnData(
            X=np.concatenate(trajectories),
            obs=self.meta.copy(),
            var=pd.DataFrame(index=population_names),
        )

    def get_coupling(self, t0, t1, covariate=None):
        """Loads a coupling for a given pair of timepoints.

        Parameters
        ----------
        t0 : int or float
            Source timepoint of the transport map.
        t1 : int of float
            Destination timepoint of the transport map.
        covariate : None or (int, int), optional
            Restrict to certain covariate values. Do not restrict if None

        Returns
        -------
        tmap : anndata.AnnData
            The transport map from t0 to t1
        """
        if t0 not in self.timepoints or t1 not in self.timepoints:
            raise ValueError("Day pair {}, {} not found".format(t0, t1))

        atomic = (
            self.day_pairs is not None and (t0, t1) in self.day_pairs
        ) or self.timepoints.index(t1) == self.timepoints.index(t0) + 1

        if not atomic and covariate is not None:
            raise ValueError("Covariate-restricted transport maps can only be atomic")

        if atomic:
            if covariate is None:
                key = (t0, t1)
            else:
                cv0, cv1 = covariate
                key = (t0, t1, str(cv0), str(cv1))
            ds_or_path = self.tmaps.get(key)
            if ds_or_path is None:
                raise ValueError("No transport map found for {}", key)
            if type(ds_or_path) is anndata.AnnData:
                return ds_or_path
            with fsspec.open(ds_or_path) as f:
                ds = anndata.read(f)
            if self.cache:
                self.tmaps[key] = ds
            return ds

        else:
            path = find_path(t0, t1, self.day_pairs, self.timepoints)
            return chain_transport_maps(self, path)

    def can_push_forward(self, *populations):
        """Checks if the populations can be pushed forward.

        Parameters
        ----------
        *populations : wot.Population
            Measure over the cells at a given timepoint to be pushed forward.

        Returns
        -------
        result : bool
            True if the populations can be pushed forward

        Raises
        ------
        ValueError
            If all populations are not in the same timepoint
        """
        return self.timepoints.index(unique_timepoint(*populations)) < len(self.timepoints) - 1

    def can_pull_back(self, *populations):
        """Checks if the populations can be pulled back.

        Parameters
        ----------
        *populations : wot.Population
            Measure over the cells at a given timepoint to be pulled back.

        Returns
        -------
        result : bool
            True if the populations can be pulled back.

        Raises
        ------
        ValueError
            If all populations are not in the same timepoint
        """
        return self.timepoints.index(unique_timepoint(*populations)) > 0

    def push_forward(self, *populations, to_time=None, normalize=True, as_list=False):
        """Pushes the population forward through the computed transport maps.

        Parameters
        ----------
        *populations : wot.Population
            Measure over the cells at a given timepoint to be pushed forward.
        to_time : int or float, optional
            Destination timepoint to push forward to.
        normalize : bool, optional, default: True
            Wether to normalize to a probability distribution or keep growth.
        as_list : bool, optional, default: False
            Wether to return a list of length 1 when a single element is passed, or a Population

        Returns
        -------
        result : wot.Population
            The push forward of the input population through the proper transport map.
            Array of populations if several populations were given as input.

        Raises
        ------
        ValueError
            If there is no further timepoint to push the population forward.
        ValueError
            If several populations are given as input but dot live in the same timepoint.

        Examples
        --------
        >>> self.push_forward(pop, to_time = 2) # -> wot.Population
        Pushing several populations at once
        >>> self.push_forward(pop1, pop2, pop3) # -> list of wot.Population
        Pulling back after pushing forward
        >>> self.pull_back(self.push_forward(pop))
        Same, but several populations at once
        >>> self.pull_back(* self.push_forward(pop1, pop2, pop3))
        """
        i = self.timepoints.index(unique_timepoint(*populations))
        j = i + 1 if to_time is None else self.timepoints.index(to_time)

        if i == -1:
            raise ValueError("Timepoint not found")
        if j == -1:
            raise ValueError("Destination timepoint not found")
        if j >= len(self.timepoints):
            raise ValueError("No further timepoints. Unable to push forward")
        if i > j:
            raise ValueError("Destination timepoint is before source. Unable to push forward")

        p = np.vstack([pop.p for pop in populations])
        while i < j:
            t0 = self.timepoints[i]
            t1 = self.timepoints[i + 1]
            tmap = self.get_coupling(t0, t1)
            p = p @ tmap.X
            if normalize:
                p = (p.T / np.sum(p, axis=1)).T
            i += 1

        result = [Population(self.timepoints[i], p[k, :]) for k in range(p.shape[0])]
        if len(result) == 1 and not as_list:
            return result[0]
        else:
            return result

    def pull_back(self, *populations, to_time=None, normalize=True, as_list=False):
        """Pulls the population back through the computed transport maps.

        Parameters
        ----------
        *populations : wot.Population
            Measure over the cells at a given timepoint to be pushed forward.
        to_time : int or float, optional
            Destination timepoint to pull back to.
        normalize : bool, optional, default: True
            Wether to normalize to a probability distribution or keep growth.
        as_list : bool, optional, default: False
            Wether to return a listof length 1 when a single element is passed, or a Population

        Returns
        -------
        result : wot.Population
            The pull back of the input population through the proper transport map.
            Array of populations if several populations were given as input.

        Raises
        ------
        ValueError
            If there is no previous timepoint to pull the population back.
        ValueError
            If several populations are given as input but dot live in the same timepoint.

        Examples
        --------
        >>> self.pull_back(pop, to_time = 0) # -> wot.Population
        Pushing several populations at once
        >>> self.pull_back(pop1, pop2, pop3) # -> list of wot.Population
        Pulling back after pushing forward
        >>> self.pull_back(self.push_forward(pop))
        Same, but several populations at once
        >>> self.pull_back(* self.push_forward(pop1, pop2, pop3))
        """
        i = self.timepoints.index(unique_timepoint(*populations))
        j = i - 1 if to_time is None else self.timepoints.index(to_time)

        if i == -1:
            raise ValueError("Timepoint not found")
        if i == 0:
            raise ValueError("No previous timepoints. Unable to pull back")
        if j == -1:
            raise ValueError("Destination timepoint not found")
        if i < j:
            raise ValueError("Destination timepoint is after source. Unable to pull back")

        p = np.vstack([pop.p for pop in populations])
        while i > j:
            t1 = self.timepoints[i]
            t0 = self.timepoints[i - 1]
            tmap = self.get_coupling(t0, t1)
            p = (tmap.X @ p.T).T
            if normalize:
                p = (p.T / np.sum(p, axis=1)).T
            i -= 1

        result = [Population(self.timepoints[i], p[k, :]) for k in range(p.shape[0])]
        if len(result) == 1 and not as_list:
            return result[0]
        else:
            return result

    def population_from_ids(self, *ids, at_time):
        """Constructs a population uniformly distributed among the ids given as input.

        Parameters
        ----------
        *ids : list of str
            The list of cell ids that belong to that population.
        at_time : int or float
            The time at which to construct the population.
            Cells that come from a different time point will be ignored.

        Returns
        -------
        *populations : list of wot.Population
            A population, uniformly distributed over the cells given as input.
            Returns None if the generated population would be empty



        Examples
        --------
        >>> cell_set = [ 'cell_1', 'cell_2', 'cell_3' ]
        >>> self.population_from_ids(cell_set) # -> wot.Population
        Multiple populations at once
        >>> multi_cell_sets = {
        >>>   'set_a': [ 'cell_a1', 'cell_a2'],
        >>>   'set_b': [ 'cell_b1', 'cell_b2'],
        >>> }
        >>> self.population_from_ids(* multi_cell_sets.values()) # -> list of wot.Population

        Notes
        -----
        The Population class is a measure over the cells at a given timepoint.
        It does not necessarily sum to 1. However, this method always returns a probability distribution over the cells of that time point.
        """
        day = float(at_time)
        df = self.meta[self.meta["day"] == day]

        def get_population(ids_el):
            cell_indices = df.index.get_indexer_for(ids_el)
            cell_indices = cell_indices[cell_indices > -1]

            if len(cell_indices) == 0:
                return None
            p = np.zeros(len(df), dtype=np.float64)
            p[cell_indices] = 1.0

            return Population(day, p)

        result = [get_population(ids_el) for ids_el in ids]
        return result


def read_transport_map_dir(transport_map_url, with_covariates=False, cache=False):
    tmaps = {}  # maps day pair to transport map url

    day_regex = "([0-9]*\\.?[0-9]+)"
    tmap_prefix = ".*"
    if with_covariates:
        pattern = re.compile(
            tmap_prefix
            + "_{}_{}_cv([a-zA-Z0-9]+)_cv([a-zA-Z0-9]+)[\\.h5ad|\\.loom|\\.zarr]".format(
                day_regex, day_regex
            )
        )
    else:
        pattern = re.compile(
            tmap_prefix + "_{}_{}[\\.h5ad|\\.loom|\\.zarr]".format(day_regex, day_regex)
        )
    scheme = get_scheme(transport_map_url)
    fs = fsspec.filesystem(scheme)

    for path in fs.glob(transport_map_url + "/*[.zarr,.h5ad]"):
        m = pattern.match(path)
        if scheme != "file":
            path = scheme + "://" + path
        if m is not None:
            try:
                t1 = float(m.group(1))
                t2 = float(m.group(2))
                if with_covariates:
                    cv1 = m.group(3)
                    cv2 = m.group(4)
                    tmaps[(t1, t2, cv1, cv2)] = path
                else:
                    tmaps[(t1, t2)] = path
            except ValueError:
                print("Unable to find day pair for " + path)
                pass

    if len(tmaps) == 0:
        raise ValueError("No transport maps found")
    day_pairs = set()
    timepoints = set()
    tmap_keys = list(tmaps.keys())
    tmap_keys.sort(key=lambda x: x[0])
    meta = None
    for i in range(len(tmap_keys)):
        key = tmap_keys[i]
        t0 = key[0]
        t1 = key[1]
        day_pairs.add((t0, t1))
        timepoints.add(t0)
        timepoints.add(t1)
        if not with_covariates:
            path = tmaps[key]
            if path.endswith(".zarr"):
                raise ValueError("Not yet implemented")
            else:
                with fs.open(path) as stream:
                    with h5py.File(stream, "r") as f:
                        # only read column ids if the last timepoint
                        if path.endswith(".loom"):
                            rids = f["/row_attrs/id"][()].astype(str)
                            cids = (
                                f["/col_attrs/id"][()].astype(str)
                                if i == len(tmap_keys) - 1
                                else None
                            )
                        else:
                            obs = f["/obs"]
                            var = f["/var"]
                            obs_key = obs.attrs.get("_index", "index")
                            var_key = var.attrs.get("_index", "index")
                            rids = obs[obs_key][()].astype(str)
                            cids = var[var_key][()].astype(str) if i == len(tmap_keys) - 1 else None
                        rdf = pd.DataFrame(index=rids, data={"day": t0})
                        cdf = (
                            pd.DataFrame(index=cids, data={"day": t1}) if cids is not None else None
                        )
                        if meta is None:
                            meta = rdf if cdf is None else pd.concat((rdf, cdf), copy=False)
                        else:
                            meta = (
                                pd.concat((meta, rdf), copy=False)
                                if cdf is None
                                else pd.concat((meta, rdf, cdf), copy=False)
                            )

    timepoints = sorted(timepoints)
    return TransportMapModel(
        tmaps=tmaps, meta=meta, timepoints=timepoints, day_pairs=day_pairs, cache=cache
    )
