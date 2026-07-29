import sys
# build shell8p0_TS_coords.xyz = TS substrate (24) + the SAME body block from shell8p0_coords.xyz
react_design=open("shell8p0_coords.xyz").read().splitlines()
n=int(react_design[0].split()[0])
body_block=react_design[2+24:2+n]   # bodies only (after the 24 substrate atoms)
ts=open("ts.xyz").read().splitlines()
ts_sub=ts[2:26]                     # TS substrate 24 atoms
out=[str(24+len(body_block)),"shell8p0 TS-under-field (TS substrate + same 8A bodies)"]+ts_sub+body_block
open("shell8p0ts_coords.xyz","w").write("\n".join(out)+"\n")
# charges identical (bodies unchanged, substrate is QM so its 24 entries are 0)
open("shell8p0ts_charges.txt","w").write(open("shell8p0_charges.txt").read())
print("wrote shell8p0ts_coords.xyz (%d atoms: 24 TS substrate + %d body) + charges"%(24+len(body_block),len(body_block)))
