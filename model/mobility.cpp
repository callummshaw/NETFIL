#include "network.h"

//radiation model for daily trips between villages
//radiaiton model from "A universal model for mobility and migration patterns" by Simini et al.
void region::radt_model(char m){
    
    //deciding which distance we want to use as basis for rad model
    double *d = NULL;
    if(m == 'r') d = road_dst;
    else if (m == 'e') d = euclid_dst;
    struct _comp_cnode_s{ //comparing distances for sorting
        bool operator() (const group::c_node *p, const group::c_node *q){ return (p->dis < q->dis);}
    } _smaller;
    
    //iterating over groups
    for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){
        group *src = j->second;

        //resetting the previous containers
        src->commuting_dist.clear();
        src->commuting_pop.clear();
        src->day_population.clear();
        src->commuting_cumsum.clear();

        src->total_commute = 0;

        if(src->commuting_dist.size() == 0){ 
            int src_id = src->gid;
            for(map<int, group*>::iterator k = groups.begin(); k != groups.end(); ++k){
                group *dst = k->second;
                int dst_id = dst->gid;

                if(dst_id == src_id) continue; //same group!

                //finding dist index!
                int index = (min(src_id, dst_id)-1)*(group_blocks*2-min(src_id, dst_id))/2 + abs(dst_id-src_id) - 1;
                src->commuting_pop.insert(pair<int, double>(dst_id, 0));
                src->commuting_dist.push_back(new group::c_node(dst_id, d[index]));                
            }
            stable_sort(src->commuting_dist.begin(), src->commuting_dist.end(), _smaller); //sorting 
        }

        double mi = src->group_pop.size(); //population of current group
        double Ti = mi*commuting_prop; //how many people will be commuting 
        double cum_sum_ceiling = 0.0;
        double com_prop;
        //now looping over all other locations
        double total_move = 0;
        double total_prop = 0;

        for(int k = 0; k < src->commuting_dist.size(); ++k){
            group *dst = groups[src->commuting_dist[k]->gid];  //other group
            double nj = dst ->group_pop.size(); //other group population
            double sij = 0; //see paper (number of people in other groups that live within radius dij (distance from current to target group), exlcuding population from i and j)
            for(int i = 0; i < k; ++i){ //iterating over other groups that have a smaller distance
                sij += groups[src->commuting_dist[i]->gid]->group_pop.size();
            } 
            com_prop = mi*nj/(mi+sij)/(mi+nj+sij);
            
            total_move += Ti*com_prop;//people from I to J ;
            total_prop += com_prop;
            src->commuting_pop[src->commuting_dist[k]->gid] = com_prop; //storing prop of people that go from i that go to j 
    
        }

        for(map<int, double>::iterator ii = src->commuting_pop.begin(); ii != src->commuting_pop.end(); ++ii){
            int gid = ii->first;
            cum_sum_ceiling += ii->second/total_prop; 
            src->commuting_cumsum[gid] = cum_sum_ceiling;
        }
        src->total_commute = total_move;
    }
}