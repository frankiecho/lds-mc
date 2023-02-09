using Distributed
addprocs(3);
@everywhere f(s,count)=(println("process id = $(myid()) s = $s count = $count");repeat(s,count))
pmap((a1,a2)->f(a1,a2),["a","b","c"],[2,1,3])