#!/usr/bin/env python3
"""Render PH-STAT's computed demonstration as a publication-quality GIF.

Requires Python 3, NumPy and Pillow. Inputs come from export_workflow_data.m.
No illustrative results are substituted for computed toolbox outputs.
"""
from pathlib import Path
import argparse, json, math, time
import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[2]
W,H,S = 1600,1000,2
BG='#F7F9FC'; WHITE='#FFFFFF'; INK='#152B45'; MUTED='#52667D'
GRID='#E1E8F0'; LIGHT='#A6B6C6'; TEAL='#008D9C'; BLUE='#246ED0'
PURPLE='#7155B6'; ORANGE='#C76D28'; PINK='#C44872'
COLORS=[TEAL,BLUE,PURPLE,ORANGE]
FONTDIR=Path('/System/Library/Fonts/Supplemental')
FONT_REG=FONTDIR/'Arial.ttf'; FONT_BOLD=FONTDIR/'Arial Bold.ttf'
if not FONT_REG.exists():
    FONT_REG=Path('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf')
    FONT_BOLD=Path('/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf')
FONTS={}

def font(size,bold=False):
    key=(size,bold)
    if key not in FONTS:
        FONTS[key]=ImageFont.truetype(str(FONT_BOLD if bold else FONT_REG),round(size*S))
    return FONTS[key]

def mix(a,b,t):
    a=np.array(tuple(bytes.fromhex(a.lstrip('#')))); b=np.array(tuple(bytes.fromhex(b.lstrip('#'))))
    return tuple(np.rint((1-t)*a+t*b).astype(int))

class Canvas:
    def __init__(self):
        self.im=Image.new('RGB',(W*S,H*S),BG)
        self.d=ImageDraw.Draw(self.im)
    def text(self,x,y,text,size=24,color=INK,bold=False,anchor='lt'):
        selected=font(size,bold)
        if any(ch in str(text) for ch in '₀₁₂') and (FONTDIR/'Arial Unicode.ttf').exists():
            key=(size,'unicode')
            if key not in FONTS:FONTS[key]=ImageFont.truetype(str(FONTDIR/'Arial Unicode.ttf'),round(size*S))
            selected=FONTS[key]
        self.d.text((round(x*S),round(y*S)),str(text),fill=color,font=selected,anchor=anchor)
    def line(self,xy,color=INK,width=2):
        self.d.line([(round(x*S),round(y*S)) for x,y in xy],fill=color,width=max(1,round(width*S)),joint='curve')
    def rect(self,box,fill=WHITE,outline=None,width=1,radius=0):
        box=tuple(round(x*S) for x in box)
        if radius:
            self.d.rounded_rectangle(box,radius*S,fill=fill,outline=outline,width=round(width*S))
        else:
            self.d.rectangle(box,fill=fill,outline=outline,width=round(width*S))
    def circle(self,x,y,r,fill=INK,outline=None,width=1):
        self.d.ellipse(tuple(round(v*S) for v in (x-r,y-r,x+r,y+r)),fill=fill,outline=outline,width=round(width*S))
    def poly(self,pts,fill):
        self.d.polygon([(round(x*S),round(y*S)) for x,y in pts],fill=fill)
    def arrow(self,a,b,color=INK,width=3,head=11):
        a=np.asarray(a,float);b=np.asarray(b,float);v=b-a
        n=np.linalg.norm(v)
        if n<1:return
        v=v/n; perp=np.array([-v[1],v[0]])
        self.line([a,b],color,width)
        self.poly([b,b-head*v+0.48*head*perp,b-head*v-0.48*head*perp],color)
    def image(self,im,box):
        x,y,w,h=box
        self.im.paste(im.resize((round(w*S),round(h*S)),Image.Resampling.NEAREST),(round(x*S),round(y*S)))
    def finish(self):
        return self.im.resize((W,H),Image.Resampling.LANCZOS)

TABS=['Morse filtration','Rips complexes','Graph filtration','Inference','Clustering','Hodge flows','Edge diffusion']
TITLES=['From a signal to persistent features','See topology emerge across scales',
        'Track network structure across thresholds','Test differences between groups',
        'Compare and organize network patterns','Separate sources, local cycles, and global flow',
        'Watch edge signals settle into harmonic circulation']
SUBS=['Smooth a noisy signal, then pair minima and maxima.',
      'Connect nearby points; fill triangles when all three edges are present.',
      'Follow connected components and cycles as weaker edges are removed.',
      'Compare topological separation with variability within groups.',
      'Represent each network by its birth and death values.',
      'Decompose an edge signal into gradient, curl, and harmonic components.',
      'Diffuse an edge signal on the same five-node simplicial complex.']
FOOT=['WFS_COS  /  PH_morse1D  /  PH_PD_display',
      'PH_rips  /  PH_boundary  /  PH_boundary_betti',
      'PH_betti  /  PH_betti_display',
      'WS_decompose  /  WS_pdist2  /  WS_ratio  /  WS_transpositions',
      'simulate_circle_difference  /  WS_decompose  /  WS_cluster  /  WS_embed',
      'Hodge_2Skeleton  /  PH_boundary  /  Hodge_decompose',
      'PH_boundary  /  Hodge_laplacian  /  spectral diffusion as in SCRIPT.mlx']

def base(scene,progress):
    c=Canvas()
    c.rect((0,0,W,9),INK)
    c.text(56,34,'PH-STAT',48,INK,True)
    c.line([(309,39),(309,94)],GRID,2)
    c.text(337,40,'Statistical inference on persistent homology',30,INK,True)
    c.text(338,78,'From signals and networks to interpretable topological analysis',21,MUTED)
    c.rect((1370,40,1544,89),'#E6F3F5',radius=7)
    c.text(1457,64,'MATLAB toolbox',20,TEAL,True,'mm')
    for i,t in enumerate(TABS):
        x=56+i*213
        c.line([(x,135),(x+194,135)],TEAL if i<scene else GRID,4)
        if i==scene:
            c.line([(x,135),(x+max(1,194*progress),135)],TEAL,4)
        c.text(x,151,f'{i+1:02d}',17,TEAL if i==scene else MUTED,True)
        c.text(x+31,151,t,17,INK if i==scene else MUTED,i==scene)
    c.text(56,207,TITLES[scene],38,INK,True)
    c.text(56,258,SUBS[scene],24,MUTED)
    c.line([(56,918),(1544,918)],GRID,2)
    c.text(56,939,FOOT[scene],19,MUTED)
    c.text(1544,939,f'{scene+1} / 7',20,TEAL,True,'rt')
    c.text(56,973,'Computed synthetic examples adapted from SCRIPT.mlx',16,MUTED,'')
    c.text(1544,973,'Moo K. Chung  ·  University of Wisconsin–Madison',16,MUTED,False,'rt')
    return c

def panel(c,box,title,kicker=None):
    c.rect(box,WHITE,GRID,1,12)
    c.text(box[0]+28,box[1]+25,title,23,INK,True)
    if kicker:c.text(box[2]-27,box[1]+29,kicker,17,MUTED,False,'rt')

def axes(c,box,xlim,ylim,xlabel,yticks,xticks,ylabel=''):
    x0,y0,x1,y1=box
    for y in yticks:
        yy=y1-(y-ylim[0])/(ylim[1]-ylim[0])*(y1-y0)
        c.line([(x0,yy),(x1,yy)],GRID,1)
        c.text(x0-13,yy,str(y),18,MUTED,False,'rm')
    for x in xticks:
        xx=x0+(x-xlim[0])/(xlim[1]-xlim[0])*(x1-x0)
        c.text(xx,y1+13,str(x),18,MUTED,False,'mt')
    c.line([(x0,y0),(x0,y1),(x1,y1)],LIGHT,1.5)
    c.text((x0+x1)/2,y1+45,xlabel,21,MUTED,False,'mt')
    if ylabel:c.text(x0,y0-30,ylabel,19,MUTED)
    def transform(points):
        p=np.asarray(points,float)
        return np.column_stack([x0+(p[:,0]-xlim[0])/(xlim[1]-xlim[0])*(x1-x0), y1-(p[:,1]-ylim[0])/(ylim[1]-ylim[0])*(y1-y0)])
    return transform

def chartline(c,tr,x,y,color=TEAL,width=3):
    pts=tr(np.column_stack([x,y]))
    if len(pts)>1:c.line(pts,color,width)

def dash(c,a,b,color=LIGHT,width=2):
    a=np.array(a);b=np.array(b);length=np.linalg.norm(b-a)
    for t in np.arange(0,1,12/max(length,1)):
        c.line([a+t*(b-a),a+min(1,t+6/max(length,1))*(b-a)],color,width)

def labeldot(c,x,y,label,color):
    c.circle(x,y,5,color);c.text(x+14,y,label,18,MUTED,False,'lm')

def caption(c,text):
    c.text(800,875,text,24,INK,False,'mm')

def sampleidx(p,n):
    return min(n-1,int(np.clip((p-0.10)/0.72,0,1)*(n-1)))

def heatmap(a,lo=None,hi=None):
    a=np.array(a,float);lo=float(np.min(a)) if lo is None else lo;hi=float(np.max(a)) if hi is None else hi
    t=np.clip((a-lo)/max(hi-lo,1e-12),0,1)
    stops=np.array([[245,248,252],[188,224,229],[42,157,169],[28,64,104]])
    f=t*3;i=np.minimum(2,f.astype(int));w=(f-i)[...,None]
    rgb=stops[i]*(1-w)+stops[i+1]*w
    return Image.fromarray(rgb.astype('uint8'))

def colorbar(c,box,lo,hi,label):
    x,y,w,h=box
    im=heatmap(np.linspace(0,1,256)[None,:],0,1)
    c.image(im,box)
    c.text(x,y+h+8,f'{lo:.2f}',16,MUTED)
    c.text(x+w,y+h+8,f'{hi:.2f}',16,MUTED,False,'rt')
    c.text(x+w/2,y+h+8,label,17,MUTED,False,'mt')

def morse(data,p):
    c=base(0,p);d=data['morse']
    x=np.array(d['x']);y=np.array(d['observed']);s=np.array(d['smooth']);pairs=np.array(d['pairs'])
    panel(c,(56,312,965,823),'Observed signal and smoothed structure')
    panel(c,(991,312,1544,823),'Persistence diagram','Finite pairs')
    tr=axes(c,(129,423,918,695),(0,1),(-0.4,2.9),'Position',[0,1,2],[0,0.5,1],'Signal amplitude')
    for xx,yy in tr(np.column_stack([x[::2],y[::2]])):
        c.circle(xx,yy,2.4,'#CFD9E3')
    n=max(2,int(len(x)*np.clip(p/0.35,0,1)))
    chartline(c,tr,x[:n],s[:n],TEAL,4)
    lim=float(np.max(pairs))+0.18
    tp=axes(c,(1072,423,1498,695),(0,lim),(0,lim),'Birth value',[0,1,2],[0,1,2],'Death value')
    dash(c,*tp(np.array([[0,0],[lim,lim]])),LIGHT,1.5)
    labeldot(c,150,794,'Noisy measurements',LIGHT)
    labeldot(c,440,794,'Smoothed signal',TEAL)
    if p>0.32:
        seq=np.argsort(pairs[:,1]);count=min(len(seq),1+int((p-.32)/.16))
        for j in seq[:count]:
            b,de=pairs[j];cx,cy=tp(np.array([[b,de]]))[0]
            c.circle(cx,cy,9,TEAL,WHITE,2)
            for val in [b,de]:
                ind=int(np.argmin(abs(s-val)));xx,yy=tr(np.array([[x[ind],s[ind]]]))[0]
                c.circle(xx,yy,7,WHITE,TEAL,2)
        c.text(1277,794,'Longer persistence = farther from diagonal',17,MUTED,False,'mm')
    caption(c,'Summarize features that persist beyond small fluctuations in the signal.')
    return c.finish()

def rips(data,p):
    c=base(1,p);d=data['rips'];idx=sampleidx(p,len(d['radii']));r=d['radii'][idx]
    points=np.array(d['points']);edges=np.array(d['edges'])-1;faces=np.array(d['faces'])-1
    eb=np.array(d['edge_birth']);fb=np.array(d['face_birth']);betti=np.array(d['betti'])
    panel(c,(56,312,909,823),'Vietoris–Rips complex',f'Radius ε = {r:.2f}')
    panel(c,(935,312,1544,823),'Topological summary')
    pts=np.column_stack([480+points[:,0]*174,565-points[:,1]*174])
    # Light fills preserve the visibility of the central hole until it closes.
    for f in faces[fb<=r]:c.poly(pts[f], '#D8EEF1')
    for a,b in edges[eb<=r]:c.line([pts[a],pts[b]],TEAL,1.6)
    for xx,yy in pts:c.circle(xx,yy,7,INK,WHITE,1.5)
    c.text(480,782,'16 points  →  edges  →  filled triangles',21,MUTED,False,'mm')
    b0,b1=betti[idx]
    c.text(983,390,str(int(b0)),54,BLUE,True)
    c.text(1060,404,'connected components  β₀',21,MUTED)
    c.text(983,470,str(int(b1)),54,PURPLE,True)
    c.text(1060,484,'unfilled loops  β₁',21,MUTED)
    tr=axes(c,(997,600,1496,724),(0,1.9),(0,1.15),'Radius ε',[0,1],[0,0.5,1,1.5],'Independent loops')
    chartline(c,tr,d['radii'][:idx+1],betti[:idx+1,1],PURPLE,4)
    xx,yy=tr([[r,b1]])[0];c.circle(xx,yy,6,PURPLE)
    message='Points begin as separate components.' if b0>1 else ('A connected ring contains one unfilled loop.' if b1 else 'Triangles fill the loop as the radius increases.')
    caption(c,message)
    return c.finish()

def graph(data,p):
    c=base(2,p);d=data['graph'];idx=sampleidx(p,len(d['thresholds']));tau=d['thresholds'][idx]
    pts0=np.array(d['coords']);C=np.array(d['weights'])
    pts=np.column_stack([430+pts0[:,0]*190,574-pts0[:,1]*165])
    panel(c,(56,312,853,823),'Weighted network',f'Keep weights > {tau:.2f}')
    panel(c,(879,312,1544,823),'Betti curves')
    for i in range(len(pts)):
        for j in range(i+1,len(pts)):
            if C[i,j]>tau:
                c.line([pts[i],pts[j]],mix('#E6EDF5',TEAL,C[i,j]),1.0+2.4*C[i,j])
    for i,(xx,yy) in enumerate(pts):
        c.circle(xx,yy,12,INK,WHITE,2)
        c.text(xx,yy,str(i+1),12,WHITE,True,'mm')
    tr0=axes(c,(956,417,1492,541),(0,1),(0,12),'',[0,6,12],[0,0.5,1],'Connected components β₀')
    tr1=axes(c,(956,633,1492,746),(0,1),(0,60),'Threshold τ',[0,30,60],[0,0.5,1],'Independent cycles β₁')
    for tr,arr,col in [(tr0,d['beta0'],BLUE),(tr1,d['beta1'],PURPLE)]:
        chartline(c,tr,d['thresholds'],arr,'#DEE5EF',2)
        chartline(c,tr,d['thresholds'][:idx+1],arr[:idx+1],col,4)
        xx,yy=tr([[tau,arr[idx]]])[0];c.circle(xx,yy,6,col)
    c.text(455,782,f'{int(d["beta0"][idx])} components     ·     {int(d["beta1"][idx])} independent cycles',22,INK,True,'mm')
    caption(c,'Increasing the threshold removes edges: components increase and cycles decrease.')
    return c.finish()

def inference(data,p):
    c=base(3,p);d=data['inference'];D=np.array(d['D01']);null=np.array(d['null']);obs=d['observed']
    panel(c,(56,312,683,823),'Between-network topological differences')
    panel(c,(709,312,1544,823),'Label-randomization test')
    c.text(364,385,'Group A: 4 networks    |    Group B: 6 networks',20,MUTED,False,'mm')
    c.image(heatmap(D,0,D.max()),(180,426,335,335))
    split=180+335*.4;h=426+335*.4
    c.line([(split,426),(split,761)],WHITE,3);c.line([(180,h),(515,h)],WHITE,3)
    c.text(247,414,'A',20,TEAL,True,'mb');c.text(414,414,'B',20,PURPLE,True,'mb')
    c.text(160,493,'A',20,TEAL,True,'rm');c.text(160,660,'B',20,PURPLE,True,'rm')
    colorbar(c,(542,426,23,335),0,D.max(),'') if False else None
    colorbar(c,(180,779,335,12),0,D.max(),'Squared W₂ distance')
    # A logarithmic x-axis keeps a large observed ratio and the null visible.
    bins=np.geomspace(max(0.1,null.min()*.85),obs*1.2,23)
    counts,_=np.histogram(null,bins=bins)
    ytop=int(np.ceil(max(counts)/10)*10)
    tr=axes(c,(788,452,1490,704),(math.log10(bins[0]),math.log10(bins[-1])),(0,ytop),'Ratio statistic R (log scale)',[0,ytop//2,ytop],[],'Number of labelings')
    k=max(1,int(len(null)*np.clip((p-.08)/.70,0,1)))
    rng=np.random.default_rng(20261016);order=rng.permutation(len(null))
    current,_=np.histogram(null[order[:k]],bins=bins)
    for i,n in enumerate(current):
        a,b=tr([[math.log10(bins[i]),0],[math.log10(bins[i+1]),n]])
        if n:c.rect((a[0]+2,b[1],b[0]-2,a[1]),TEAL)
    for v in [.5,1,2,5,10,20,50]:
        if bins[0]<=v<=bins[-1]:
            xx=tr([[math.log10(v),0]])[0,0];c.text(xx,716,str(v),18,MUTED,False,'mt')
    xx=tr([[math.log10(obs),0]])[0,0]
    dash(c,(xx,449),(xx,704),PINK,2.5)
    c.text(1485,385,f'Observed R = {obs:.2f}',22,PINK,True,'rt')
    c.text(790,385,f'{k} / {len(null)} label allocations',20,MUTED)
    if k==len(null):
        c.text(1125,793,f'Exact upper-tail p = {d["p_exact"]:.4f}',25,INK,True,'mm')
    else:c.text(1125,793,'Build the reference distribution by relabeling groups',21,MUTED,False,'mm')
    caption(c,'Synthetic demonstration: all 210 allocations of 10 networks into groups of 4 and 6.')
    return c.finish()

def clustering(data,p):
    c=base(4,p);d=data['clustering'];xy=np.array(d['xy']);groups=np.array(d['groups'])
    panel(c,(56,312,717,823),'Four simulated network groups')
    panel(c,(743,312,1544,823),'Topological embedding','32 networks')
    for j,cloud in enumerate(d['clouds']):
        cloud=np.array(cloud);cx=229+(j%2)*316;cy=489+(j//2)*220
        lo=cloud.min(0);hi=cloud.max(0);norm=(cloud-(lo+hi)/2)/max(hi-lo)*195
        for q in norm:c.circle(cx+q[0],cy-q[1],3.2,COLORS[j])
        c.text(cx,cy+103,f'Group {j+1}',18,COLORS[j],True,'mm')
    mins=xy.min(0);maxs=xy.max(0);pad=(maxs-mins)*.15
    xt=np.round(np.linspace(mins[0],maxs[0],3),2);yt=np.round(np.linspace(mins[1],maxs[1],3),2)
    tr=axes(c,(844,439,1486,700),(mins[0]-pad[0],maxs[0]+pad[0]),(mins[1]-pad[1],maxs[1]+pad[1]),'Mean birth value',yt,xt,'Mean death value')
    n=max(1,int(len(xy)*np.clip((p-.08)/.62,0,1)))
    order=np.array([i+8*j for i in range(8) for j in range(4)])
    for q in order[:n]:
        x,y=tr(xy[q:q+1])[0];c.circle(x,y,7,COLORS[groups[q]-1],WHITE,1.5)
    c.text(1143,793,'Colors indicate known generating groups',20,MUTED,False,'mm')
    caption(c,'Cluster using full birth/death features; view networks through a two-dimensional summary.')
    return c.finish()

def flowgraph(c,data,box,values,color,phase,scale=None,labels=True):
    x0,y0,x1,y1=box;coords=np.array(data['coords']);edges=np.array(data['edges'])-1
    pts=np.column_stack([(x0+x1)/2+coords[:,0]*(x1-x0)*.43, (y0+y1)/2-coords[:,1]*(y1-y0)*.40])
    faces=np.array(data['faces'],int).reshape(-1,3)-1
    for face in faces:c.poly(pts[face],'#EAF0F5')
    scale=max(abs(np.array(values))) if scale is None else scale
    for i,(a,b) in enumerate(edges):
        v=float(values[i]);c.line([pts[a],pts[b]],GRID,2)
        if abs(v)<1e-7:continue
        start,end=(pts[a],pts[b]) if v>0 else (pts[b],pts[a]);delta=end-start;unit=delta/np.linalg.norm(delta)
        start=start+16*unit;end=end-18*unit
        weight=abs(v)/max(scale,1e-10)
        width=1.8+4.2*weight
        c.arrow(start,end,color,width,11)
        # Markers are direction cues, not particle simulations.
        for off in [0,.5]:
            t=(phase*.65+off+i*.13)%1;pos=start+t*(end-start)
            c.circle(*pos,3+1.2*weight,WHITE,color,1.2)
    for i,(x,y) in enumerate(pts):
        c.circle(x,y,13,WHITE,INK,2)
        if labels:c.text(x,y,str(i+1),13,INK,True,'mm')

def hodge(data,p):
    c=base(5,p);d=data['hodge'];names=['Input flow','Gradient','Curl','Harmonic']
    keys=['input','gradient','curl','harmonic'];cols=[INK,BLUE,ORANGE,PURPLE]
    desc=['Observed edge signal','Source–sink variation','Filled-triangle circulation','Unfilled-cycle circulation']
    scale=max(abs(np.array(d['input'])))
    for j in range(4):
        x=56+j*380
        panel(c,(x,354,x+348,791),names[j])
        c.rect((x+1,354,x+347,360),cols[j])
        flowgraph(c,d,(x+23,443,x+325,696),d[keys[j]],cols[j],p*8,scale)
        c.text(x+174,751,desc[j],18,MUTED,False,'mm')
    for x,sign in [(418,'='),(798,'+'),(1178,'+')]:c.text(x,563,sign,31,MUTED,True,'mm')
    c.text(800,325,'X = XG + XC + XH',24,INK,True,'mm')
    caption(c,'Arrow width shows flow magnitude. Moving markers indicate flow direction.')
    return c.finish()

def diffusion(data,p):
    c=base(6,p);d=data['hodge'];idx=sampleidx(p,len(d['time']));t=d['time'][idx]
    flows=np.array(d['diffusion']);tt=np.array(d['time']);Xh=np.array(d['harmonic'])
    panel(c,(56,312,751,823),'Evolving edge flow',f'Diffusion time t = {t:.2f}')
    panel(c,(777,312,1544,823),'Transient signals decay; harmonic flow remains')
    flowgraph(c,d,(158,441,650,735),flows[:,idx],PURPLE,p*10,max(abs(np.array(d['input']))))
    c.text(403,788,'Filled triangle is included in the Hodge Laplacian',18,MUTED,False,'mm')
    tr=axes(c,(865,441,1490,700),(0,5),(-.65,2.15),'Diffusion time t',[0,1,2],[0,1,2,3,4,5],'Signed edge flow')
    cols=[TEAL,BLUE,PURPLE,ORANGE,PINK,'#5D7787']
    for i,y in enumerate(flows):
        chartline(c,tr,tt,y,'#E6ECF2',1.5)
        chartline(c,tr,tt[:idx+1],y[:idx+1],cols[i],3)
        xx,yy=tr([[t,y[idx]]])[0];c.circle(xx,yy,4,cols[i])
    c.text(1177,792,'Each curve represents one edge',20,MUTED,False,'mm')
    caption(c,'Hodge diffusion removes gradient and curl components and retains harmonic circulation.')
    return c.finish()

def overview(data):
    c=Canvas()
    c.rect((0,0,W,10),INK)
    c.text(56,48,'PH-STAT',64,INK,True)
    c.text(1544,67,'MATLAB  /  BRAIN NETWORK ANALYSIS',20,TEAL,True,'rt')
    c.text(56,135,'Statistical inference on persistent homology',36,INK,True)
    c.text(56,193,'From signals and network structure to statistical evidence',29,MUTED)
    names=['Signals & networks','Topological structure','Network comparisons','Network flows']
    details=[('Signals and connectivity','Morse filtration · smoothing'),
             ('Components and cycles','Rips complexes · Betti curves'),
             ('Distances and group patterns','Inference · clustering · embedding'),
             ('Sources and circulation','Hodge decomposition · diffusion')]
    for i,name in enumerate(names):
        x=56+380*i
        c.rect((x,285,x+348,761),WHITE,GRID,1,12)
        c.rect((x+1,286,x+347,292),COLORS[i])
        c.text(x+24,321,name,24,INK,True)
        c.text(x+174,686,details[i][0],19,INK,True,'mm')
        c.text(x+174,724,details[i][1],16,MUTED,False,'mm')
    C=np.array(data['graph']['weights'])
    c.image(heatmap(C,0,1),(96,400,268,242))
    c.text(230,373,'Weighted connectivity',18,MUTED,False,'mm')
    rd=data['rips'];pts=np.array(rd['points']);pts=np.c_[610+pts[:,0]*120,517-pts[:,1]*120]
    for face,t in zip(np.array(rd['faces'])-1,rd['face_birth']):
        if t<.85:c.poly(pts[face],'#DDEEF6')
    for edge,t in zip(np.array(rd['edges'])-1,rd['edge_birth']):
        if t<.85:c.line(pts[edge],BLUE,2)
    for x,y in pts:c.circle(x,y,5.5,INK,WHITE,1)
    c.text(610,373,'A loop across scales',18,MUTED,False,'mm')
    cd=data['clustering'];xy=np.array(cd['xy']);lo=xy.min(0);hi=xy.max(0)
    xy=(xy-lo)/(hi-lo)
    c.line([(865,405),(865,642),(1118,642)],LIGHT,1.5)
    for q,g in zip(xy,cd['groups']):c.circle(883+q[0]*222,625-q[1]*210,5.4,COLORS[g-1],WHITE,1)
    c.text(990,373,'Topological embedding',18,MUTED,False,'mm')
    hd=data['hodge'];flowgraph(c,hd,(1230,423,1505,657),hd['harmonic'],PURPLE,.9,None)
    c.text(1370,373,'Harmonic circulation',18,MUTED,False,'mm')
    c.text(56,813,'Explore the worked examples in SCRIPT.mlx',31,INK,True)
    c.text(56,861,'Filtrations  →  topological summaries  →  statistical analysis',25,TEAL)
    c.line([(56,920),(1544,920)],GRID,2)
    c.text(56,947,'Computed synthetic examples  ·  Seven animated modules',20,MUTED)
    c.text(1544,947,'Moo K. Chung  ·  University of Wisconsin–Madison',19,MUTED,False,'rt')
    return c.finish()

SCENES=[morse,rips,graph,inference,clustering,hodge,diffusion]

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--preview',action='store_true')
    args=parser.parse_args()
    data_path=ROOT/'tools/workflow/data.json'
    if not data_path.exists():data_path=ROOT/'build/workflow/data.json'
    data=json.loads(data_path.read_text())
    out=ROOT/'assets/workflow';out.mkdir(parents=True,exist_ok=True)
    build=ROOT/'build/workflow';build.mkdir(parents=True,exist_ok=True)
    keyframes=[]
    for i,render in enumerate(SCENES):
        im=render(data,.88 if i!=6 else .93)
        im.save(build/f'scene-{i+1}.png')
        keyframes.append(im)
    sheet=Image.new('RGB',(1600,4*520),BG)
    for i,im in enumerate(keyframes):
        thumb=im.resize((800,500),Image.Resampling.LANCZOS)
        sheet.paste(thumb,((i%2)*800,(i//2)*520))
    sheet.save(build/'contact-sheet.png')
    poster=overview(data)
    poster.save(out/'ph-stat-workflow-overview.png')
    if args.preview:
        print('Preview images saved',flush=True);return
    # One palette for all frames: consistent colors and efficient GIF deltas.
    atlas=Image.new('RGB',(400*4,250*2),BG)
    for i,im in enumerate(keyframes+[poster]):atlas.paste(im.resize((400,250)),((i%4)*400,(i//4)*250))
    palette=atlas.quantize(colors=192,method=Image.Quantize.MEDIANCUT,dither=Image.Dither.NONE)
    palette_rgb=palette.getpalette()[:192*3]
    for col in [INK,MUTED,LIGHT,TEAL,BLUE,PURPLE,ORANGE,PINK,BG,WHITE,GRID,'#5D7787']:
        for tint in [0,.25,.5,.75]:palette_rgb.extend(mix(col,WHITE,tint))
    palette_rgb.extend([255]*(768-len(palette_rgb)))
    palette.putpalette(palette_rgb)
    frames=[poster.quantize(palette=palette,dither=Image.Dither.NONE)];durations=[2500]
    for i,render in enumerate(SCENES):
        print(f'Rendering {i+1}/7: {TABS[i]}',flush=True)
        # 120 frames at 20 fps, plus a short final hold for each module.
        for j in range(120):
            im=render(data,j/119)
            frames.append(im.quantize(palette=palette,dither=Image.Dither.NONE))
            durations.append(50)
        durations[-1]=650
    gif=out/'ph-stat-workflow.gif'
    frames[0].save(gif,save_all=True,append_images=frames[1:],duration=durations,
        loop=0,optimize=True,disposal=1)
    with Image.open(gif) as chk:
        count=chk.n_frames;duration=0
        for i in range(count):
            chk.seek(i);duration+=chk.info.get('duration',0)
    info={'width':W,'height':H,'frames':count,'duration_seconds':duration/1000,
          'nominal_fps':20,'bytes':gif.stat().st_size,'megabytes_decimal':round(gif.stat().st_size/1e6,3),
          'size_limit_bytes':25000000,'loop':0,'matlab_source_commit':data['source_commit']}
    assert gif.stat().st_size<25000000,info
    (build/'gif-validation.json').write_text(json.dumps(info,indent=2)+'\n')
    print(json.dumps(info,indent=2),flush=True)

if __name__=='__main__':main()
