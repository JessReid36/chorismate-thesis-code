import sys
# usage: make_ts_any.py <tag e.g. shell5p0>  (needs <tag>_coords.xyz, <tag>_charges.txt, ts.xyz in dir)
tag=sys.argv[1]
design=open("%s_coords.xyz"%tag).read().splitlines()
n=int(design[0].split()[0])
body_block=design[2+24:2+n]
ts=open("ts.xyz").read().splitlines(); ts_sub=ts[2:26]
out=[str(24+len(body_block)),"%s TS-under-field"%tag]+ts_sub+body_block
open("%sts_coords.xyz"%tag,"w").write("\n".join(out)+"\n")
open("%sts_charges.txt"%tag,"w").write(open("%s_charges.txt"%tag).read())
print("wrote %sts_coords.xyz (%d atoms) + charges"%(tag,24+len(body_block)))
