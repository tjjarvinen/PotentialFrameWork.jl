module PotentialFrameWorkNeighbourListsExt
import NeighbourLists
using PotentialFrameWork.NeighbourList: NeighbourSnapshot, neighbours

# Dispatch only on our type; upstream constructors/types remain upstream-owned.
NeighbourLists.neighbours(list::NeighbourSnapshot,i::Integer)=neighbours(list,i)
NeighbourLists.nsites(list::NeighbourSnapshot)=length(list.indices)
NeighbourLists.npairs(list::NeighbourSnapshot)=sum(length,list.indices;init=0)
end
