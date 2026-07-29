import sys
# build shell7p4_TS_coords.xyz = TS substrate (24) + the SAME body block from shell7p4_coords.xyz
react_design=open("shell7p4_coords.xyz").read().splitlines()
n=int(react_design[0].split()[0])
body_block=react_design[2+24:2+n]   # bodies only (after the 24 substrate atoms)
ts=open("ts.xyz").read().splitlines()
ts_sub=ts[2:26]                     # TS substrate 24 atoms
out=[str(24+len(body_block)),"shell7p4 TS-under-field (TS substrate + same 8A bodies)"]+ts_sub+body_block
open("shell7p4ts_coords.xyz","w").write("\n".join(out)+"\n")
# charges identical (bodies unchanged, substrate is QM so its 24 entries are 0)
open("shell7p4ts_charges.txt","w").write(open("shell7p4_charges.txt").read())
print("wrote shell7p4ts_coords.xyz (%d atoms: 24 TS substrate + %d body) + charges"%(24+len(body_block),len(body_block)))
