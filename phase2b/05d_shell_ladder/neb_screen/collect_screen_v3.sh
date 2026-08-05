#!/bin/bash
# GOCAT gate, CORRECT per-endpoint integrity:
#  R valid = converged + ether intact (O3-C4<2.2)
#  P valid = converged + C1-C6 formed (<1.9) + not shattered (O3-C4<5.0; ref product O3-C4=2.9)
echo "design | R:O3C4 C1C6 conv ok | P:O3C4 C1C6 conv ok | verdict"
echo "-------|--------------------|--------------------|-------"
for tag in 2p0 7p0 7p1 7p2 7p3 7p4 7p5 7p6 7p7 7p8 7p9 8p0 8p1 8p2 8p3 8p5 8p6 8p7 8p8 8p9 9p0 10p0 11p0 12p0 13p0 15p0; do
  rx=run_screen_${tag}_reactant/shell${tag}_reactant_screened.xyz
  px=run_screen_${tag}_product/shell${tag}_product_screened.xyz
  rlog=run_screen_${tag}_reactant/scr_${tag}_reactant.log
  plog=run_screen_${tag}_product/scr_${tag}_product.log
  g(){ python3 -c "
import math
try:
 L=open('$1').read().splitlines();a=[[float(v) for v in L[2+i].split()[1:4]] for i in range(24)]
 d=lambda i,j: math.sqrt(sum((a[i][k]-a[j][k])**2 for k in range(3)))
 print('%.3f %.3f'%(d(7,8),d(0,12)))
except: print('--- ---')
"; }
  c(){ [ -f "$1" ] && { grep -qi "optimization converged" "$1" && echo 1 || { grep -qi "maximum number" "$1" && echo 0 || echo R; }; } || echo -; }
  [ -f "$rx" ] || [ -f "$px" ] || { printf '%-6s| running\n' "$tag"; continue; }
  rg=$(g "$rx"); pg=$(g "$px"); rc=$(c "$rlog"); pc=$(c "$plog")
  python3 -c "
rg='$rg'.split(); pg='$pg'.split(); rc='$rc'; pc='$pc'
def rok():
    try: return rc=='1' and float(rg[0])<2.2
    except: return False
def pok():
    try: return pc=='1' and float(pg[1])<1.9 and float(pg[0])<5.0
    except: return False
R=rok(); P=pok()
rflag='Y' if R else 'N'; pflag='Y' if P else 'N'
if rc in('R','-') or pc in('R','-'): v='running/partial'
elif R and P: v='PASS -> NEB'
else:
    f=[]
    if rc=='0': f.append('R-noconv')
    try:
        if float(rg[0])>=2.2: f.append('R-etherbroken')
    except: pass
    if pc=='0': f.append('P-noconv')
    try:
        if float(pg[1])>=1.9: f.append('P-noC1C6')
        if float(pg[0])>=5.0: f.append('P-shattered')
    except: pass
    v='FAIL('+','.join(f)+')' if f else 'FAIL'
print('%-6s| %s %s   %s  %s | %s %s   %s  %s | %s'%('$tag',rg[0],rg[1],rc,rflag,pg[0],pg[1],pc,pflag,v))
"
done
