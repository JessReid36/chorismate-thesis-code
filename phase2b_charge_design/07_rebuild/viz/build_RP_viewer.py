#!/usr/bin/env python3
"""Standalone HTML viewer of reactant & product with atom index (0-based) + element labels."""
import os, json
HERE = os.path.dirname(os.path.abspath(__file__))
INP  = os.path.join(HERE, "..", "inputs")

def load(name):
    L = open(os.path.join(INP, name)).read().splitlines()
    n = int(L[0].split()[0]); out = []
    for ln in L[2:2+n]:
        p = ln.split(); out.append([p[0], float(p[1]), float(p[2]), float(p[3])])
    return out

geo = {"reactant": load("reactant.xyz"), "product": load("product.xyz")}

TEMPLATE = r'''<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>Chorismate R / P - atom numbering</title>
<style>
 body{margin:0;font-family:system-ui,sans-serif;background:#0e0e12;color:#ddd;overflow:hidden}
 #bar{position:fixed;top:0;left:0;right:0;padding:8px 12px;background:#16161c;z-index:10;
      display:flex;gap:14px;align-items:center;box-shadow:0 1px 6px #0008;flex-wrap:wrap}
 button{background:#26262f;color:#ddd;border:1px solid #3a3a46;border-radius:6px;padding:6px 12px;cursor:pointer}
 button.active{background:#3a5;border-color:#3a5;color:#fff}
 label{font-size:13px;user-select:none}.hint{font-size:12px;color:#888}
 canvas{display:block}#reacting{color:#ff6}
</style></head><body>
<div id="bar">
 <button id="btnR" class="active">Reactant</button>
 <button id="btnP">Product</button>
 <label><input type="checkbox" id="idx" checked> atom #</label>
 <label><input type="checkbox" id="elem" checked> element</label>
 <label><input type="checkbox" id="bonds" checked> bonds</label>
 <span class="hint">drag=rotate | scroll=zoom | <span id="reacting">O3-C4 (7-8) &amp; C1-C6 (0-12) yellow</span></span>
</div>
<canvas id="c"></canvas>
<script>
const GEO=__GEO__;
const COLORS={C:'#909090',H:'#e8e8e8',O:'#ff4030',N:'#3050f0'};
const RCOV={C:0.77,H:0.31,O:0.66,N:0.70};
const REACTING=new Set(['7-8','0-12']);
let cur='reactant',rotX=0.4,rotY=0.6,zoom=1,dragging=false,lx=0,ly=0;
const cv=document.getElementById('c'),ctx=cv.getContext('2d');
function resize(){cv.width=innerWidth;cv.height=innerHeight;}resize();onresize=()=>{resize();draw();};
function center(a){let c=[0,0,0];a.forEach(x=>{c[0]+=x[1];c[1]+=x[2];c[2]+=x[3]});c=c.map(v=>v/a.length);
 return a.map(x=>[x[0],x[1]-c[0],x[2]-c[1],x[3]-c[2]]);}
function rot(p){let[x,y,z]=p;let cyv=Math.cos(rotY),syv=Math.sin(rotY);[x,z]=[x*cyv-z*syv,x*syv+z*cyv];
 let cxv=Math.cos(rotX),sxv=Math.sin(rotX);[y,z]=[y*cxv-z*sxv,y*sxv+z*cxv];return[x,y,z];}
function draw(){
 ctx.clearRect(0,0,cv.width,cv.height);
 const atoms=center(GEO[cur]);
 const showIdx=document.getElementById('idx').checked,showEl=document.getElementById('elem').checked,showB=document.getElementById('bonds').checked;
 const S=48*zoom,cx=cv.width/2,cy=cv.height/2+20;
 const proj=atoms.map(a=>{const r=rot([a[1],a[2],a[3]]);return{el:a[0],x:cx+r[0]*S,y:cy-r[1]*S,z:r[2]};});
 if(showB){for(let i=0;i<atoms.length;i++)for(let j=i+1;j<atoms.length;j++){
   const a=atoms[i],b=atoms[j],d=Math.hypot(a[1]-b[1],a[2]-b[2],a[3]-b[3]),cut=1.3*((RCOV[a[0]]||.7)+(RCOV[b[0]]||.7));
   if(d<cut){const key=i+'-'+j;ctx.strokeStyle=REACTING.has(key)?'#ff6':'#555';ctx.lineWidth=REACTING.has(key)?3:1.5;
    ctx.beginPath();ctx.moveTo(proj[i].x,proj[i].y);ctx.lineTo(proj[j].x,proj[j].y);ctx.stroke();}}}
 const order=[...proj.keys()].sort((p,q)=>proj[p].z-proj[q].z);
 for(const i of order){const p=proj[i],rad=(p.el==='H'?7:11)*Math.sqrt(zoom);
  ctx.beginPath();ctx.arc(p.x,p.y,rad,0,7);ctx.fillStyle=COLORS[p.el]||'#c0f';ctx.fill();
  ctx.lineWidth=1;ctx.strokeStyle='#000';ctx.stroke();
  if(showIdx||showEl){ctx.fillStyle='#fff';ctx.font='bold 12px system-ui';ctx.textAlign='left';
   let lbl=(showEl?p.el:'')+(showIdx?i:'');ctx.fillText(lbl,p.x+rad+1,p.y-2);}}
 ctx.fillStyle='#888';ctx.font='12px system-ui';ctx.textAlign='left';
 ctx.fillText(cur.toUpperCase()+'  ('+atoms.length+' atoms) - label = element+index (0-based)',12,cv.height-14);
}
document.getElementById('btnR').onclick=()=>{cur='reactant';btnR.classList.add('active');btnP.classList.remove('active');draw();};
document.getElementById('btnP').onclick=()=>{cur='product';btnP.classList.add('active');btnR.classList.remove('active');draw();};
['idx','elem','bonds'].forEach(id=>document.getElementById(id).onchange=draw);
cv.onmousedown=e=>{dragging=true;lx=e.clientX;ly=e.clientY;};
onmouseup=()=>dragging=false;
onmousemove=e=>{if(dragging){rotY+=(e.clientX-lx)*0.01;rotX+=(e.clientY-ly)*0.01;lx=e.clientX;ly=e.clientY;draw();}};
cv.onwheel=e=>{e.preventDefault();zoom*=e.deltaY<0?1.1:0.9;zoom=Math.max(0.3,Math.min(6,zoom));draw();};
draw();
</script></body></html>'''
open(os.path.join(HERE,"view_RP.html"),"w").write(TEMPLATE.replace("__GEO__", json.dumps(geo)))
print("wrote view_RP.html")
for name in ("reactant","product"):
    print("\n== %s ==" % name)
    for i,a in enumerate(geo[name]):
        tag = "  <-- reacting" if i in (0,7,8,12) else ""
        print("  %2d  %-2s  %8.3f %8.3f %8.3f%s" % (i,a[0],a[1],a[2],a[3],tag))
