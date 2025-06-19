//Meth
function sin(x)     {return Math.sin(x);}
function cos(x)     {return Math.cos(x);}
function cosh(x)    {return Math.cosh(x);}
function tan(x)     {return Math.tan(x);}
function tanh(x)    {return Math.tanh(x);}

function qsign(x) { return typeof x === 'number' ? x ? x < 0 ? -1 : 1 : x === x ? 0 : NaN : NaN; }

function qvec3minx(vec3a, vec3b) { return Math.min(vec3a[0], vec3b[0]); }
function qvec3miny(vec3a, vec3b) { return Math.min(vec3a[1], vec3b[1]); }
function qvec3minz(vec3a, vec3b) { return Math.min(vec3a[2], vec3b[2]); }

function qvec3maxx(vec3a, vec3b) { return Math.max(vec3a[0], vec3b[0]); }
function qvec3maxy(vec3a, vec3b) { return Math.max(vec3a[1], vec3b[1]); }
function qvec3maxz(vec3a, vec3b) { return Math.max(vec3a[2], vec3b[2]); }

function qtrunc_signed(a) { return a|0; }
function norm(min, max, x) { return (x - min) / (max - min); }

function matMul(mat1, mat2)
{
    const m1lx = mat1.length
    const m1ly = mat1[0].length
    const m2lx = mat2.length
    const m2ly = mat2[0].length
    
    //Invalid matrix check
    if 
    (
        m1lx <= 0 || m1ly <= 0 || m2lx <= 0 || m2ly <= 0 ||
        m1ly != m2lx
    
    ) return null;
    
    //Prepare return matrix
    let ret = new Array(m1lx);
    for (let i = 0; i < m1lx; ++i) 
    {
        ret[i] = new Array(m2ly).fill(0);
    }

    //Do multiplication
    /*
                            [11b,12b,13b]
        [11a,12a,13a]  x    [21b,22b,23b]   =   [(11a*11b + 12a*21b + 13a * 31b),   (11a*12b + 12a*22b + 13a * 32b),    (11a*13b + 12a*23b + 13a * 33b)]
                            [31b,32b,33b]
    */
    
    
    for (let y = 0; y < m1lx; ++y) 
    {
        for (let x = 0; x < m2ly; ++x) 
        {
            for (let k = 0; k < m1ly; ++k) 
            {
                ret[y][x] += mat1[y][k] * mat2[k][x];  
            }
        }
    }

    return ret;
}

function lerpsl(sl, start, end)
{
    const x1 = sl[start];
    const x2 = sl[end];
    const t = 1/(start - end);

    for (let i = start; i < end; ++i) 
    {
        sl[i] = (x2 - x1) * t + x1;
    }
}

//fast linear interpolation
//Return vec3
function lerpvec3_f(ax, ay, az, bx, by, bz, t)
{
    return [
                ((bx-ax)*t+ax),
                ((by-ay)*t+ay),
                ((bz-az)*t+az)
            ];
}

//Buffers

//vertex:
/*
    [ vert_id, vec2_size, vec3_offset, float[y][x] depth]
*/
function make_vbuff(id, sx, sy, offx, offy, offz)
{
    let d1 = new Array(sy);
    for (let i = 0; i < sy; ++i) 
    {
        d1[i] = new Array(sx).fill(NaN);
    }
    return [
        id,
        [sx, sy],
        [offx, offy, offz],
        d1
    ];
} 

//Rasterization
function qrast_line(vec3a, vec3b, mp, xpad, ypad, zpad, out_vbuff, quality) 
{
    for (let i = 0.0; i <= 1.0; i += 1/quality) 
    {
        const k = lerpvec3_f(vec3a[0], vec3a[1], vec3a[2], vec3b[0], vec3b[1], vec3b[2], i);
        const x = qtrunc_signed(mp*k[0]+mp+xpad);
        const y = qtrunc_signed(mp*k[1]+mp+ypad);
        const z = qtrunc_signed(mp*k[2]+mp+zpad);
        if (x < 0 || y < 0 || x >= out_vbuff[1][0] || y >= out_vbuff[1][1]) continue;
        
        out_vbuff[3][y][x] = z; 
    }
}

//Compositing
function composite_tri_sides_equ(vecbuff_s_1, vecbuff_s_2, vecbuff_s_3)
{
    const sx = vecbuff_s_1[1][0];
    const sy = vecbuff_s_1[1][1];

    for (let y = 0; y < sy; ++y) 
    {    
        for (let x = 0; x < sx; ++x) 
        {        
            vecbuff_s_1[3][y][x] = (!Number.isNaN(vecbuff_s_1[3][y][x]) && vecbuff_s_1[3][y][x] > vecbuff_s_2[3][y][x]) ? vecbuff_s_1[3][y][x] : (!Number.isNaN(vecbuff_s_2[3][y][x]) && vecbuff_s_2[3][y][x] > vecbuff_s_3[3][y][x]) ? vecbuff_s_2[3][y][x] : vecbuff_s_3[3][y][x];  
        }
    }   
}


function qrast_tri(id, vec3a, vec3b, vec3c, w, h, _ratio, _xpad, _ypad, quality) 
{
    if (!vec3a || !vec3b || !vec3c)
    {
        return null;
    }

    //minmax
    const minx = Math.min(vec3a[0], vec3b[0], vec3c[0]);
    const miny = Math.min(vec3a[1], vec3b[1], vec3c[1]);
    const minz = Math.min(vec3a[2], vec3b[2], vec3c[2]);

    const maxx = Math.max(vec3a[0], vec3b[0], vec3c[0]);
    const maxy = Math.max(vec3a[1], vec3b[1], vec3c[1]);
    const maxz = Math.max(vec3a[2], vec3b[2], vec3c[2]);

    //diff
    const diffx = qtrunc_signed((maxx - minx)*w);
    const diffy = qtrunc_signed((maxy - miny)*h);
    const diffz = (maxy - miny);
    const _mp = _ratio / 2;

    //Everything is in NSPC baka!

    let b1 = make_vbuff(id, diffx, diffy, minx, miny, minz);

    for (let y = 0; y < diffy; ++y) 
    {
        for (let x = 0; x < diffx; ++x) 
        {
            b1[3][y][x] = 0;
        }
    }
    
    return b1;
}

//2d rendering
function render2DLine_noDepth(dat, canv_w, canv_h, mp, xpad, ypad, vec3a, vec3b, rgbacolor, quality) 
{
    /*
        -1.0 -> mp*-1 + mp = 0
        0.0 -> mp* 0 + mp = mp
        1.0 -> mp* 1 + mp = 2mp = ratio
        
        apply padding so stays in middle
        z clip at < 0
    */

    for (let i = 0.0; i <= 1.0; i += 1/quality) 
    {
        const k = lerpvec3_f(vec3a[0], vec3a[1], vec3a[2], vec3b[0], vec3b[1], vec3b[2], i);
        //if (qtrunc_signed(k[2]) < 0) continue;
        const x = qtrunc_signed(mp*k[0]+mp+xpad);
        const y = qtrunc_signed(mp*k[1]+mp+ypad);
        if (x < 0 || y < 0 || x >= canv_w || y >= canv_h) continue;
        const pos = 4*(y*canv_w+x);
        

        dat.data[pos]        = rgbacolor[0];
        dat.data[pos + 1]    = rgbacolor[1];
        dat.data[pos + 2]    = rgbacolor[2];
        dat.data[pos + 3]    = rgbacolor[3];
    }
    return dat;
}

function render2DTriFace_noDepth(dat, canv_w, canv_h, mp, xpad, ypad, vec3a, vec3b, vec3c, rgbacolor) 
{
    for (let i = 0.0; i <= 1.0; i += 1/quality) 
    {
        const k = lerpvec3_f(vec3a[0], vec3a[1], vec3a[2], vec3b[0], vec3b[1], vec3b[2], i);
        //if (qtrunc_signed(k[2]) < 0) continue;
        const x = qtrunc_signed(mp*k[0]+mp+xpad);
        const y = qtrunc_signed(mp*k[1]+mp+ypad);
        if (x < 0 || y < 0 || x >= canv_w || y >= canv_h) continue;
        const pos = 4*(y*canv_w+x);
        

        dat.data[pos]        = rgbacolor[0];
        dat.data[pos + 1]    = rgbacolor[1];
        dat.data[pos + 2]    = rgbacolor[2];
        dat.data[pos + 3]    = rgbacolor[3];
    }
    return dat;
}

//3d rendering
function renderTri(w, h, vecbuff_tri, rgbacolor, out_buff) 
{
    for (let y = 0; y < vecbuff_tri[1][1]; ++y) 
    {
        for (let x = 0; x < vecbuff_tri[1][0]; ++x) 
        {
            if (x < 0 || y < 0 || x >= w || y >= h) continue;
            if (Number.isNaN(vecbuff_tri[3][y][x])) continue;
            const pos = 4*(y*w+x);

            out_buff.data[pos]        = rgbacolor[0];
            out_buff.data[pos + 1]    = rgbacolor[1];
            out_buff.data[pos + 2]    = rgbacolor[2];
            out_buff.data[pos + 3]    = rgbacolor[3];
        }   
    }
}

//Matrix constatnts
const proj_mat = 
[
    [   1,  0,  0   ],
    [   0,  1,  0   ],
    [   0,  0,  1   ]
];

//θ
const θ = 0.01;
const rotX = 
[
    [   1,  0,  0   ],
    [   0, Math.cos(θ), -Math.sin(θ)    ],
    [   0, Math.sin(θ),  Math.cos(θ)    ]
];

const rotY =
[
    [   Math.cos(θ), 0, Math.sin(θ) ],
    [   0,  1,  0  ],
    [  -Math.sin(θ), 0, Math.cos(θ) ]
];

const rotZ =
[
    [   Math.cos(θ), -Math.sin(θ), 0    ],
    [   Math.sin(θ), Math.cos(θ), 0     ],
    [   0,  0,  1   ]
];


function drawCtx(ctx, dat) {ctx.putImageData(dat, 0, 0);}

const canv = document.getElementById("canv");

canv.width = window.innerWidth;
canv.height = window.innerHeight;

var _ctx = canv.getContext("2d", { willReadFrequently: true });
if (_ctx == null) 
{
    _ctx = canv.getContext("2d");
}
const ctx = _ctx;

const w = canv.width;
const h = canv.height;

const ratio = (w < h) ? h:w;
const padX = (w-ratio) / 2;
const padY = (h-ratio) / 2;

const midpoint = ratio / 2;

// -1.0 - 1.0
//Clip mask : drop z < 0
let pointsTri = 
[
    // x      y       z
    [0.0,    0.0,   0.0],
    [-1.0,   0.0,    0.0],
    [1.0,   -1.0,   0.0]
];

const color_1 = [255, 0, 0, 255];
const color_2 = [0, 255, 0, 255];
const color_3 = [0, 0, 255, 255];
const quality = 200;

let lastTime = 0;
const targetFPS = 120;
const frameInterval = 1000 / targetFPS;

function render(currentTime)
{
    const deltaTime = currentTime - lastTime;
    
    if (deltaTime > frameInterval)
    {
        ctx.fillStyle = "black";
        ctx.fillRect(0, 0, w, h);
        const imageData = ctx.getImageData(0, 0, w, h);

        for (let i = 0; i < pointsTri.length; i++)
        {
            pointsTri[i] = matMul([pointsTri[i]], rotX)[0];
            pointsTri[i] = matMul([pointsTri[i]], rotY)[0];
            pointsTri[i] = matMul([pointsTri[i]], rotZ)[0];
        }
        
        render2DLine_noDepth(imageData, w, h, midpoint, padX, padY, pointsTri[0], pointsTri[1], color_1, quality);
        render2DLine_noDepth(imageData, w, h, midpoint, padX, padY, pointsTri[1], pointsTri[2], color_2, quality);
        render2DLine_noDepth(imageData, w, h, midpoint, padX, padY, pointsTri[2], pointsTri[0], color_3, quality);
        
        ctx.putImageData(imageData, 0, 0);
        //lastTime = currentTime - (deltaTime % frameInterval);
    }

    requestAnimationFrame(render);
    
}

let vecbuff;
function renderTriFilled(currentTime)
{
    const deltaTime = currentTime - lastTime;
    
    if (deltaTime > frameInterval)
    {
        ctx.fillStyle = "black";
        ctx.fillRect(0, 0, w, h);
        const imageData = ctx.getImageData(0, 0, w, h);
        
        for (let i = 0; i < pointsTri.length; i++)
        {
            pointsTri[i] = matMul([pointsTri[i]], rotX)[0];
            pointsTri[i] = matMul([pointsTri[i]], rotY)[0];
            pointsTri[i] = matMul([pointsTri[i]], rotZ)[0];
        }

        vecbuff = qrast_tri(0, pointsTri[0], pointsTri[1], pointsTri[2], w, h, ratio, padX, padY, 200);
        renderTri(w, h, vecbuff, color_1, imageData);

        ctx.putImageData(imageData, 0, 0);
        render(currentTime);
        lastTime = currentTime - (deltaTime % frameInterval);
    }
    
    requestAnimationFrame(renderTriFilled);
}   
requestAnimationFrame(render);