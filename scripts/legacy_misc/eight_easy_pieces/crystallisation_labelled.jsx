import { useState, useEffect, useRef, useCallback } from "react";

const project = (x, y, z, cx, cy, scale, rotY, rotX) => {
  const cosY = Math.cos(rotY), sinY = Math.sin(rotY);
  const cosX = Math.cos(rotX), sinX = Math.sin(rotX);
  let x1 = x * cosY - z * sinY;
  let z1 = x * sinY + z * cosY;
  let y1 = y * cosX - z1 * sinX;
  let z2 = y * sinX + z1 * cosX;
  const p = 600 / (600 + z2);
  return { px: cx + x1 * scale * p, py: cy + y1 * scale * p, depth: z2 };
};

// Q3 edges
const Q3_EDGES = [];
for (let i = 0; i < 8; i++)
  for (let j = i + 1; j < 8; j++) {
    const xor = i ^ j;
    if (xor === 1 || xor === 2 || xor === 4) Q3_EDGES.push([i, j]);
  }

const OCTA_V = [
  [1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]
];
const OCTA_EDGES = [];
for (let a = 0; a < 6; a++)
  for (let b = a + 1; b < 6; b++)
    if (!((a===0&&b===1)||(a===2&&b===3)||(a===4&&b===5)))
      OCTA_EDGES.push([a, b]);
const OCTA_FACES = [
  [0,2,4],[0,2,5],[0,3,4],[0,3,5],
  [1,2,4],[1,2,5],[1,3,4],[1,3,5]
];

// Face centroids
const FACE_C = [];
for (let i = 0; i < 8; i++) {
  FACE_C.push([
    ((i & 4) ? -1 : 1) / 3,
    ((i & 2) ? -1 : 1) / 3,
    ((i & 1) ? -1 : 1) / 3
  ]);
}

// Bit labels and their antipodal pair colours
const BIT_INFO = [
  { label: "G\u2080", pair: 0 },   // 000 - pair with W (111)
  { label: "G\u2081", pair: 1 },   // 001 - pair with chi (110)
  { label: "LQ",      pair: 2 },   // 010 - pair with I3 (101)
  { label: "C\u2080", pair: 3 },   // 011 - pair with C1 (100)
  { label: "C\u2081", pair: 3 },   // 100 - pair with C0 (011)
  { label: "I\u2083", pair: 2 },   // 101 - pair with LQ (010)
  { label: "\u03C7",  pair: 1 },   // 110 - pair with G1 (001)
  { label: "W",       pair: 0 },   // 111 - pair with G0 (000)
];

// Antipodal pair colours (matching the paper figures)
const PAIR_COLORS = [
  { r: 100, g: 140, b: 210, name: "Blue" },     // G0 <-> W
  { r: 165, g: 115, b: 195, name: "Purple" },    // G1 <-> chi
  { r: 70,  g: 185, b: 155, name: "Teal" },      // LQ <-> I3
  { r: 210, g: 90,  b: 70,  name: "Red/Green" }, // C0 <-> C1
];

const CLUSTER_COLORS = [
  {r:231,g:76,b:60},{r:46,g:204,b:113},{r:52,g:152,b:219},
  {r:243,g:156,b:18},{r:155,g:89,b:182},{r:26,g:188,b:156},
  {r:230,g:126,b:34},{r:200,g:200,b:210},
];

const PHASES = [
  "Random qubit soup \u2014 no structure",
  "Entanglement bonds breaking and reforming\u2026",
  "Clusters nucleating from energy minimisation\u2026",
  "Parity-check circuits closing on Q\u2083 faces\u2026",
  "Octahedral voids crystallising\u2026",
  "Voids aligning \u2014 bridge edges forming\u2026",
  "Lattice crystallisation complete"
];

function getLatticePositions(n) {
  const sp = 4.0;
  if (n === 1) return [[0,0,0]];
  if (n === 2) return [[-sp/2,0,0],[sp/2,0,0]];
  if (n === 3) return [[-sp*0.5,0,-sp*0.29],[sp*0.5,0,-sp*0.29],[0,0,sp*0.58]];
  if (n === 4) return [[-sp/2,0,-sp/2],[sp/2,0,-sp/2],[-sp/2,0,sp/2],[sp/2,0,sp/2]];
  const pos = [];
  const r = sp * n / (2 * Math.PI) * 1.1;
  for (let i = 0; i < n; i++) {
    const a = (2*Math.PI*i)/n;
    pos.push([Math.cos(a)*r, 0, Math.sin(a)*r]);
  }
  return pos;
}

function findBridges(cps, s) {
  const bridges = [];
  for (let ci = 0; ci < cps.length; ci++)
    for (let cj = ci+1; cj < cps.length; cj++) {
      let best = null, bestD = 5.5;
      for (let vi = 0; vi < 6; vi++)
        for (let vj = 0; vj < 6; vj++) {
          const dx = (cps[ci][0]+OCTA_V[vi][0]*s)-(cps[cj][0]+OCTA_V[vj][0]*s);
          const dy = (cps[ci][1]+OCTA_V[vi][1]*s)-(cps[cj][1]+OCTA_V[vj][1]*s);
          const dz = (cps[ci][2]+OCTA_V[vi][2]*s)-(cps[cj][2]+OCTA_V[vj][2]*s);
          const d = Math.sqrt(dx*dx+dy*dy+dz*dz);
          if (d < bestD) { bestD = d; best = {ci,cj,vi,vj,d}; }
        }
      if (best) bridges.push(best);
    }
  return bridges;
}

export default function App() {
  const canvasRef = useRef(null);
  const animRef = useRef(null);
  const stateRef = useRef(null);
  const [phase, setPhase] = useState(0);
  const [isRunning, setIsRunning] = useState(false);
  const [nQubits, setNQubits] = useState(24);
  const [speed, setSpeed] = useState(1);
  const [stats, setStats] = useState(null);

  const initState = useCallback((n) => {
    const nodes = [];
    for (let i = 0; i < n; i++)
      nodes.push({
        x:(Math.random()-0.5)*7, y:(Math.random()-0.5)*7, z:(Math.random()-0.5)*7,
        cluster:-1, faceIndex:-1, targetX:0, targetY:0, targetZ:0,
      });
    const edges=[], edgeSet=new Set(), deg=new Array(n).fill(0);
    for (let i = 0; i < n; i++) {
      const c=[];
      for (let j=0;j<n;j++) {
        if(i!==j&&deg[j]<3){const k=`${Math.min(i,j)}-${Math.max(i,j)}`;if(!edgeSet.has(k))c.push(j);}
      }
      while(deg[i]<3&&c.length>0){
        const idx=Math.floor(Math.random()*c.length), j=c[idx];
        if(deg[j]<3){const k=`${Math.min(i,j)}-${Math.max(i,j)}`;if(!edgeSet.has(k)){edges.push([i,j]);edgeSet.add(k);deg[i]++;deg[j]++;}}
        c.splice(idx,1);
      }
    }
    return {nodes,edges,clusters:[],frame:0,maxFrames:550,edgesReplaced:false,bridgesComputed:false,bridges:[],clusterPositions:[],alignPhase:false};
  },[]);

  const startAnim = useCallback(() => {
    stateRef.current = initState(nQubits);
    setPhase(0); setStats(null); setIsRunning(true);
  }, [nQubits, initState]);

  useEffect(() => {
    if (!isRunning) return;
    const canvas = canvasRef.current;
    const ctx = canvas.getContext("2d");
    const W=canvas.width, H=canvas.height, cx=W/2, cy=H/2;
    let rotY = 0;
    const octS = 1.6;

    const draw = () => {
      const S = stateRef.current;
      if (!S) return;
      S.frame += speed;
      rotY += 0.005 * speed;
      const rotX = 0.30;
      const t = S.frame / S.maxFrames;
      const numCl = Math.floor(nQubits/8);
      const leftover = nQubits - numCl*8;

      let cp;
      if(t<0.10) cp=0; else if(t<0.22) cp=1; else if(t<0.40) cp=2;
      else if(t<0.55) cp=3; else if(t<0.72) cp=4;
      else if(t<0.90) cp=5; else cp=6;
      setPhase(cp);

      // Assign clusters
      if(t>0.22 && S.clusters.length===0){
        const rem=[...Array(S.nodes.length).keys()];
        for(let i=rem.length-1;i>0;i--){const j=Math.floor(Math.random()*(i+1));[rem[i],rem[j]]=[rem[j],rem[i]];}
        const clusters=[];
        const sR=numCl<=3?3.2:2.8, aS=(2*Math.PI)/Math.max(numCl,1);
        for(let c=0;c<numCl;c++){
          const cl=[], ox=Math.cos(aS*c)*sR, oz=Math.sin(aS*c)*sR;
          for(let i=0;i<8;i++){const ni=rem.shift();cl.push(ni);S.nodes[ni].cluster=c;S.nodes[ni].faceIndex=i;
            const fc=FACE_C[i];S.nodes[ni].targetX=ox+fc[0]*octS;S.nodes[ni].targetY=fc[1]*octS;S.nodes[ni].targetZ=oz+fc[2]*octS;}
          clusters.push(cl);
        }
        for(const idx of rem) S.nodes[idx].cluster=-2;
        S.clusters=clusters;
        S.clusterPositions=clusters.map((_,ci)=>[Math.cos(aS*ci)*sR,0,Math.sin(aS*ci)*sR]);
      }

      // Replace edges
      if(t>0.45&&!S.edgesReplaced){
        S.edgesReplaced=true;
        const ne=[];
        S.clusters.forEach(cl=>{Q3_EDGES.forEach(([a,b])=>{if(cl[a]!==undefined&&cl[b]!==undefined)ne.push([cl[a],cl[b]]);});});
        S.edges=ne;
      }

      // Align
      if(t>0.70&&!S.alignPhase){
        S.alignPhase=true;
        const lp=getLatticePositions(numCl);
        S.clusterPositions=lp;
        S.clusters.forEach((cl,ci)=>{const p=lp[ci];cl.forEach((ni,fi)=>{const fc=FACE_C[fi];
          S.nodes[ni].targetX=p[0]+fc[0]*octS;S.nodes[ni].targetY=p[1]+fc[1]*octS;S.nodes[ni].targetZ=p[2]+fc[2]*octS;});});
      }

      // Bridges
      if(t>0.80&&!S.bridgesComputed&&S.clusterPositions.length>1){
        S.bridgesComputed=true;
        S.bridges=findBridges(S.clusterPositions,octS);
      }

      // Animate
      const settling=t>0.26, sStr=Math.min(1,(t-0.26)/0.3);
      const aStr=t>0.70?Math.min(1,(t-0.70)/0.15):0;
      const eR=settling?(0.03+aStr*0.04)*speed:0;
      S.nodes.forEach((node,idx)=>{
        if(settling&&node.cluster>=0){node.x+=(node.targetX-node.x)*eR;node.y+=(node.targetY-node.y)*eR;node.z+=(node.targetZ-node.z)*eR;}
        else if(node.cluster===-2){node.x+=(Math.random()-0.5)*0.06*speed;node.y+=(Math.random()-0.5)*0.06*speed;node.z+=(Math.random()-0.5)*0.06*speed;}
        else{node.x+=(Math.random()-0.5)*0.05*speed;node.y+=(Math.random()-0.5)*0.05*speed;node.z+=(Math.random()-0.5)*0.05*speed;}
      });

      // === DRAW ===
      ctx.fillStyle="#06060c";
      ctx.fillRect(0,0,W,H);

      // Q3 edges (very dim)
      S.edges.forEach(([i,j])=>{
        const ni=S.nodes[i],nj=S.nodes[j];
        const pi=project(ni.x,ni.y,ni.z,cx,cy,65,rotY,rotX);
        const pj=project(nj.x,nj.y,nj.z,cx,cy,65,rotY,rotX);
        let alpha=0.08;
        if(ni.cluster>=0&&ni.cluster===nj.cluster&&settling) alpha=Math.min(0.12,sStr*0.12);
        ctx.strokeStyle=`rgba(120,120,140,${alpha})`;
        ctx.lineWidth=0.5;
        ctx.beginPath();ctx.moveTo(pi.px,pi.py);ctx.lineTo(pj.px,pj.py);ctx.stroke();
      });

      // Octahedron wireframes (bright)
      if(t>0.50){
        const wA=Math.min(0.9,(t-0.50)*2.2);
        const fA=Math.min(0.08,(t-0.60)*0.25);
        S.clusters.forEach((cl,ci)=>{
          const col=CLUSTER_COLORS[ci%CLUSTER_COLORS.length];
          let ax=0,ay=0,az=0;
          cl.forEach(ni=>{ax+=S.nodes[ni].x;ay+=S.nodes[ni].y;az+=S.nodes[ni].z;});
          ax/=8;ay/=8;az/=8;
          const pV=OCTA_V.map(v=>project(ax+v[0]*octS,ay+v[1]*octS,az+v[2]*octS,cx,cy,65,rotY,rotX));

          // Face fills
          if(fA>0.01) OCTA_FACES.forEach(face=>{
            const pts=face.map(vi=>pV[vi]);
            ctx.beginPath();ctx.moveTo(pts[0].px,pts[0].py);ctx.lineTo(pts[1].px,pts[1].py);ctx.lineTo(pts[2].px,pts[2].py);ctx.closePath();
            ctx.fillStyle=`rgba(${col.r},${col.g},${col.b},${fA})`;ctx.fill();
          });

          // Edges (bright, thick)
          OCTA_EDGES.forEach(([a,b])=>{
            ctx.strokeStyle=`rgba(${col.r},${col.g},${col.b},${wA})`;
            ctx.lineWidth=2.0;
            ctx.beginPath();ctx.moveTo(pV[a].px,pV[a].py);ctx.lineTo(pV[b].px,pV[b].py);ctx.stroke();
          });

          // Vertices (bright dots)
          pV.forEach(pv=>{
            ctx.beginPath();ctx.arc(pv.px,pv.py,3,0,Math.PI*2);
            ctx.fillStyle=`rgba(${col.r},${col.g},${col.b},${wA})`;ctx.fill();
            ctx.beginPath();ctx.arc(pv.px-1,pv.py-1,1.2,0,Math.PI*2);
            ctx.fillStyle=`rgba(255,255,255,${wA*0.4})`;ctx.fill();
          });
        });
      }

      // Bridges (golden, brightest)
      if(S.bridges.length>0&&t>0.82){
        const bA=Math.min(1,(t-0.82)*3.5);
        S.bridges.forEach(br=>{
          let ai={x:0,y:0,z:0},aj={x:0,y:0,z:0};
          S.clusters[br.ci].forEach(ni=>{ai.x+=S.nodes[ni].x;ai.y+=S.nodes[ni].y;ai.z+=S.nodes[ni].z;});
          ai.x/=8;ai.y/=8;ai.z/=8;
          S.clusters[br.cj].forEach(ni=>{aj.x+=S.nodes[ni].x;aj.y+=S.nodes[ni].y;aj.z+=S.nodes[ni].z;});
          aj.x/=8;aj.y/=8;aj.z/=8;

          const v1=OCTA_V[br.vi],v2=OCTA_V[br.vj];
          const p1=project(ai.x+v1[0]*octS,ai.y+v1[1]*octS,ai.z+v1[2]*octS,cx,cy,65,rotY,rotX);
          const p2=project(aj.x+v2[0]*octS,aj.y+v2[1]*octS,aj.z+v2[2]*octS,cx,cy,65,rotY,rotX);

          // Glow
          ctx.shadowColor=`rgba(255,215,80,${bA*0.6})`;ctx.shadowBlur=10;
          ctx.strokeStyle=`rgba(255,225,120,${bA})`;ctx.lineWidth=3;
          ctx.beginPath();ctx.moveTo(p1.px,p1.py);ctx.lineTo(p2.px,p2.py);ctx.stroke();
          ctx.shadowBlur=0;

          // Bridge vertex highlights
          [p1,p2].forEach(p=>{
            ctx.beginPath();ctx.arc(p.px,p.py,4,0,Math.PI*2);
            ctx.fillStyle=`rgba(255,235,150,${bA*0.8})`;ctx.fill();
          });

          // Midpoint label
          const mx=(p1.px+p2.px)/2, my=(p1.py+p2.py)/2;
          if(bA>0.5){
            ctx.fillStyle=`rgba(255,225,120,${(bA-0.3)*1.2})`;
            ctx.font="bold 10px monospace";ctx.textAlign="center";
            ctx.fillText("bridge",mx,my-10);
          }
        });
      }

      // Qubit nodes with bit identity labels
      const pN=S.nodes.map((n,i)=>({...project(n.x,n.y,n.z,cx,cy,65,rotY,rotX),idx:i,node:n}));
      pN.sort((a,b)=>a.depth-b.depth);

      pN.forEach(({px,py,idx,node})=>{
        let r=60,g=60,b=75,radius=2.8;
        if(node.cluster>=0&&settling){
          const col=CLUSTER_COLORS[node.cluster%CLUSTER_COLORS.length];
          r=col.r;g=col.g;b=col.b;radius=4;
        } else if(node.cluster===-2){
          r=255;g=50;b=50;radius=2.5+Math.sin(S.frame*0.12+idx)*1;
        }

        // Glow
        const grad=ctx.createRadialGradient(px,py,0,px,py,radius*2.5);
        grad.addColorStop(0,`rgba(${r},${g},${b},0.2)`);
        grad.addColorStop(1,`rgba(${r},${g},${b},0)`);
        ctx.beginPath();ctx.arc(px,py,radius*2.5,0,Math.PI*2);ctx.fillStyle=grad;ctx.fill();

        // Body
        ctx.beginPath();ctx.arc(px,py,radius,0,Math.PI*2);
        ctx.fillStyle=`rgb(${r},${g},${b})`;ctx.fill();

        // BIT IDENTITY LABELS (always visible once crystallised)
        if(node.cluster>=0 && node.faceIndex>=0 && t>0.62){
          const labelAlpha = Math.min(0.95, (t-0.62)*3);
          const bi = BIT_INFO[node.faceIndex];
          const pc = PAIR_COLORS[bi.pair];

          // Label background pill
          const lbl = bi.label;
          const tw = ctx.measureText ? 0 : 0; // approximate
          ctx.font = "bold 9px monospace";
          const metrics = ctx.measureText(lbl);
          const pillW = metrics.width + 6;
          const pillH = 12;
          const lx = px - pillW/2;
          const ly = py - radius - pillH - 3;

          ctx.fillStyle = `rgba(${pc.r},${pc.g},${pc.b},${labelAlpha*0.25})`;
          ctx.beginPath();
          ctx.roundRect(lx, ly, pillW, pillH, 3);
          ctx.fill();

          ctx.fillStyle = `rgba(${pc.r},${pc.g},${pc.b},${labelAlpha})`;
          ctx.font = "bold 9px monospace";
          ctx.textAlign = "center";
          ctx.textBaseline = "middle";
          ctx.fillText(lbl, px, ly + pillH/2);
          ctx.textBaseline = "alphabetic";
        }
      });

      // === HUD ===
      ctx.fillStyle="rgba(6,6,12,0.85)";
      ctx.fillRect(0,H-54,W,54);
      ctx.fillStyle="#b0adc0";
      ctx.font="bold 13px system-ui, -apple-system, sans-serif";
      ctx.textAlign="center";
      ctx.fillText(PHASES[cp],cx,H-30);

      const bW=W-50;
      ctx.fillStyle="rgba(255,255,255,0.06)";ctx.fillRect(25,H-11,bW,3);
      ctx.fillStyle=cp>=6?"#2ecc71":cp>=5?"#f0c040":"#3498db";
      ctx.fillRect(25,H-11,bW*Math.min(1,t),3);

      // Info panel
      ctx.fillStyle="#555";ctx.font="10px monospace";ctx.textAlign="left";
      const energy=cp>=6?Math.round(-612*(numCl/3)):Math.round(-60-552*(numCl/3)*Math.min(1,t*1.1));
      ctx.fillText(`N=${nQubits}  E=${energy}`,10,18);
      if(S.clusters.length>0){
        ctx.fillText(`${S.clusters.length}\u00D78 octahedra`,10,32);
        if(S.bridges.length>0&&t>0.82){ctx.fillStyle="#e0c840";ctx.fillText(`${S.bridges.length} bridge${S.bridges.length!==1?"s":""}`,10,46);}
        if(leftover>0){ctx.fillStyle="#e55";ctx.fillText(`+${leftover} frustrated`,10,S.bridges.length?60:46);}
      }

      // Right info
      if(t>0.55){
        const iA=Math.min(0.4,(t-0.55));
        ctx.fillStyle=`rgba(255,255,255,${iA})`;ctx.font="9px monospace";ctx.textAlign="right";
        ctx.fillText("8 qubits on 8 triangular faces",W-10,18);
        ctx.fillText("Q\u2083 face-adjacency (dim lines)",W-10,30);
        ctx.fillText("Octahedron skeleton (bright lines)",W-10,42);
        if(S.bridges.length>0&&t>0.84){
          ctx.fillStyle=`rgba(255,225,120,${Math.min(0.5,(t-0.84)*2)})`;
          ctx.fillText("Bridge edges (golden) = gauge links",W-10,56);
        }
      }

      // Antipodal pair legend (bottom right, after crystallisation)
      if(t>0.75){
        const lgA=Math.min(0.7,(t-0.75)*2);
        const lgX=W-10, lgY=H-100;
        ctx.font="8px monospace";ctx.textAlign="right";
        ctx.fillStyle=`rgba(200,200,210,${lgA*0.6})`;
        ctx.fillText("Antipodal pairs:",lgX,lgY);
        const pairLabels = [
          "G\u2080 \u2194 W",
          "G\u2081 \u2194 \u03C7",
          "LQ \u2194 I\u2083",
          "C\u2080 \u2194 C\u2081",
        ];
        pairLabels.forEach((pl,i)=>{
          const pc=PAIR_COLORS[i];
          ctx.fillStyle=`rgba(${pc.r},${pc.g},${pc.b},${lgA})`;
          ctx.fillText(pl,lgX,lgY+13+i*12);
        });
      }

      if(cp>=6&&!stats) setStats({clusters:numCl,leftover,bridges:S.bridges.length,perfect:leftover===0});

      if(S.frame<S.maxFrames*1.15) animRef.current=requestAnimationFrame(draw);
      else setIsRunning(false);
    };

    animRef.current=requestAnimationFrame(draw);
    return ()=>{if(animRef.current)cancelAnimationFrame(animRef.current);};
  },[isRunning,nQubits,speed,stats]);

  return (
    <div className="flex flex-col items-center bg-gray-950 min-h-screen p-3">
      <h1 className="text-lg font-bold text-gray-200 mb-0.5">Quantum Vacuum Crystallisation</h1>
      <p className="text-xs text-gray-500 mb-2">Qubits self-organise into octahedral voids with labelled bit identities</p>

      <canvas ref={canvasRef} width={720} height={470} className="rounded-lg border border-gray-800 mb-3"/>

      <div className="flex flex-wrap gap-1.5 items-center justify-center mb-2">
        <span className="text-gray-500 text-xs">N:</span>
        {[16,24,32,23,25].map(n=>(
          <button key={n} onClick={()=>{setNQubits(n);setIsRunning(false);setStats(null);}}
            className={`px-2 py-1 rounded text-xs font-mono ${nQubits===n?"bg-blue-600 text-white":"bg-gray-800 text-gray-400 hover:bg-gray-700"}`}
          >{n}{n%8!==0?"*":""}</button>
        ))}
        <span className="text-gray-600 text-xs ml-2">Speed:</span>
        {[0.5,1,2].map(s=>(
          <button key={s} onClick={()=>setSpeed(s)}
            className={`px-2 py-1 rounded text-xs ${speed===s?"bg-blue-600 text-white":"bg-gray-800 text-gray-400 hover:bg-gray-700"}`}
          >{s}\u00D7</button>
        ))}
      </div>

      <button onClick={startAnim}
        className="px-6 py-2 bg-blue-600 hover:bg-blue-500 text-white rounded-lg font-medium text-sm"
      >{isRunning?"\u21BB Restart":"\u25B6 Crystallise"}</button>

      {stats && (
        <div className={`mt-3 p-3 rounded-lg border text-center max-w-lg ${stats.perfect?"bg-green-950/50 border-green-800":"bg-amber-950/50 border-amber-700"}`}>
          <div className={`text-base font-bold ${stats.perfect?"text-green-400":"text-amber-400"}`}>
            {stats.perfect?"\u2713 Lattice Crystallised":"\u26A0 Frustrated Vacuum"}
          </div>
          <div className="text-gray-400 text-xs mt-1">
            {stats.clusters} octahedr{stats.clusters===1?"on":"a"}
            {" \u2022 "}{stats.clusters*8} face-qubits (G\u2080 G\u2081 LQ C\u2080 C\u2081 I\u2083 \u03C7 W)
            {stats.bridges>0&&<span className="text-yellow-400">{" \u2022 "}{stats.bridges} bridge{stats.bridges!==1?"s":""}</span>}
            {stats.leftover>0&&<span className="text-red-400">{" \u2022 "}{stats.leftover} frustrated</span>}
          </div>
        </div>
      )}

      <div className="flex flex-wrap gap-4 mt-3 text-xs text-gray-500 justify-center">
        <div className="flex items-center gap-1.5">
          <div className="w-4 h-px bg-gray-600 opacity-30"></div><span>Q\u2083 edges (dim)</span>
        </div>
        <div className="flex items-center gap-1.5">
          <div className="w-4 h-0.5 bg-red-400 opacity-80"></div><span>Octahedron (bright)</span>
        </div>
        <div className="flex items-center gap-1.5">
          <div className="w-4 h-0.5 bg-yellow-300"></div><span>Bridge (golden)</span>
        </div>
        {PAIR_COLORS.map((pc,i)=>(
          <div key={i} className="flex items-center gap-1">
            <div className="w-2.5 h-2.5 rounded-sm" style={{backgroundColor:`rgb(${pc.r},${pc.g},${pc.b})`}}></div>
            <span>{["G\u2080\u2194W","G\u2081\u2194\u03C7","LQ\u2194I\u2083","C\u2080\u2194C\u2081"][i]}</span>
          </div>
        ))}
      </div>
    </div>
  );
}
