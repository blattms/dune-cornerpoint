/*
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services.
  Copyright 2015 Statoil AS

  This file is part of The Open Porous Media project  (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <opm/grid/common/ZoltanPartition.hpp>

#include <opm/grid/utility/OpmParserIncludes.hpp>
#include <sstream>

#if defined(HAVE_ZOLTAN) && defined(HAVE_MPI)
namespace Dune
{
namespace cpgrid
{
std::pair<std::vector<int>, std::unordered_set<std::string> >
zoltanGraphPartitionGridOnRoot(const CpGrid& cpgrid,
                               const std::vector<const OpmWellType*> * wells,
                               const double* transmissibilities,
                               const CollectiveCommunication<MPI_Comm>& cc,
                               int root)
{
    int rc = ZOLTAN_OK - 1;
    float ver = 0;
    struct Zoltan_Struct *zz;
    int changes, numGidEntries, numLidEntries, numImport, numExport;
    ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
    int *importProcs, *importToPart, *exportProcs, *exportToPart;
    int argc=0;
    char** argv = 0 ;
    bool doLoadBalance = true;
    std::vector<int> parts(cpgrid.numCells(), cc.rank());
    std::vector<std::vector<int> > wells_on_proc;

    rc = Zoltan_Initialize(argc, argv, &ver);

    if ( rc == ZOLTAN_FATAL )
    {
        OPM_THROW(std::runtime_error, "Zoltan initialization failed with fatal error.");
    }
    if ( rc == ZOLTAN_MEMERR )
    {
        OPM_THROW(std::runtime_error, "Zoltan initialization failed because of error in memory allocation.");
    }

    // In principal it should be safe to use all processors during loadbalancing
    // But for some OpenMPI version errors have been encountered if the whole
    // graph is stored on one process and the rest has no vertices and no edges.
    // Therefore we will use only on process for loadbalancing
#ifdef LOADBALANCE_GRAPH_ALL_PROCS
    zz = Zoltan_Create(cc);
#else
    if ( cc.rank() == 0 )
    {
        zz = Zoltan_Create(MPI_COMM_SELF);
    }
    else
    {
        doLoadBalance = false;
    }
#endif

    if ( doLoadBalance )
    {
        Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
        Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");
        Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
        Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");
        Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
        Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");
#ifndef LOADBALANCE_GRAPH_ALL_PROCS
	std::ostringstream partitions;
        partitions << cc.size();
        auto no = partitions.str().c_str();
        Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS", no); // Number of partition requested. Needed for MPI_COMM_SELF
        Zoltan_Set_Param(zz, "NUM_LOCAL_PARTS", no); // Number of partition requested. Needed for MPI_COMM_SELF
#endif
#ifndef NDEBUG
        Zoltan_Set_Param(zz, "CHECK_GRAPH", "2");
#else
        Zoltan_Set_Param(zz, "CHECK_GRAPH", "0");
#endif
        Zoltan_Set_Param(zz,"EDGE_WEIGHT_DIM","0");
        Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");
        Zoltan_Set_Param(zz, "PHG_EDGE_SIZE_THRESHOLD", ".35");  /* 0-remove all, 1-remove none */

        // For the load balancer one process has the whole grid and
        // all others an empty partition before loadbalancing.
        bool partitionIsEmpty     = cc.rank()!=root;
        bool partitionIsWholeGrid = !partitionIsEmpty;

        std::shared_ptr<CombinedGridWellGraph> grid_and_wells;

        if( wells )
        {
            Zoltan_Set_Param(zz,"EDGE_WEIGHT_DIM","1");
            grid_and_wells.reset(new CombinedGridWellGraph(cpgrid,
                                                           wells,
                                                           transmissibilities,
                                                           partitionIsEmpty));
            Dune::cpgrid::setCpGridZoltanGraphFunctions(zz, *grid_and_wells,
                                                        partitionIsEmpty);
        }
        else
        {
            Dune::cpgrid::setCpGridZoltanGraphFunctions(zz, cpgrid, partitionIsEmpty);
        }

        rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
                                 &changes,        /* 1 if partitioning was changed, 0 otherwise */
                                 &numGidEntries,  /* Number of integers used for a global ID */
                                 &numLidEntries,  /* Number of integers used for a local ID */
                                 &numImport,      /* Number of vertices to be sent to me */
                                 &importGlobalGids,  /* Global IDs of vertices to be sent to me */
                                 &importLocalGids,   /* Local IDs of vertices to be sent to me */
                                 &importProcs,    /* Process rank for source of each incoming vertex */
                                 &importToPart,   /* New partition for each incoming vertex */
                                 &numExport,      /* Number of vertices I must send to other processes*/
                                 &exportGlobalGids,  /* Global IDs of the vertices I must send */
                                 &exportLocalGids,   /* Local IDs of the vertices I must send */
                                 &exportProcs,    /* Process to which I send each of the vertices */
                                 &exportToPart);  /* Partition to which each vertex will belong */
        if ( rc == ZOLTAN_FATAL )
        {
            OPM_THROW(std::runtime_error, "Loadbalancing failed with fatal error.");
        }
        if ( rc == ZOLTAN_MEMERR )
        {
            OPM_THROW(std::runtime_error, "Loadbalancing failed because of error in memory allocation.");
        }

        for ( int i=0; i < numExport; ++i )
        {
 #ifdef LOADBALANCE_GRAPH_ALL_PROCS
            parts[exportLocalGids[i]] = exportProcs[i];
#else
            parts[exportLocalGids[i]] = exportToPart[i];
#endif
        }

        if( wells && partitionIsWholeGrid )
        {
            wells_on_proc =
                postProcessPartitioningForWells(parts,
                                                *wells,
                                                grid_and_wells->getWellConnections(),
                                                cc.size());

#ifndef NDEBUG
            int index = 0;
            for( auto well : grid_and_wells->getWellsGraph() )
            {
                int part=parts[index];
                std::set<std::pair<int,int> > cells_on_other;
                for( auto vertex : well )
                {
                    if( part != parts[vertex] )
                    {
                        cells_on_other.insert(std::make_pair(vertex, parts[vertex]));
                    }
                }
                if ( cells_on_other.size() )
                {
                    OPM_THROW(std::domain_error, "Well is distributed between processes, which should not be the case!");
                }
                ++index;
            }
#endif
        }
        // free space allocated for zoltan.
        Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, &exportProcs, &exportToPart);
        Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, &importProcs, &importToPart);
        Zoltan_Destroy(&zz);
    }

    std::unordered_set<std::string> defunct_well_names;

    if( wells )
    {
        defunct_well_names = computeDefunctWellNames(wells_on_proc,
                                                     *wells,
                                                     cc,
                                                     root);
    }

    cc.broadcast(&parts[0], parts.size(), root);

    return std::make_pair(parts, defunct_well_names);
}
}
}
#endif // HAVE_ZOLTAN
