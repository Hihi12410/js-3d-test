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

function qtrunc_signed(a) { return a|0; }

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

function render2DLine_noBlend(dat, canv_w, canv_h, mp, xpad, ypad, vec3a, vec3b, rgbacolor, quality) 
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

const ratio = (w > h) ? h:w;
const padX = (w-ratio) / 2;
const padY = (h-ratio) / 2;

const midpoint = ratio / 2;

// -1.0 - 1.0
//Clip mask : drop z < 0
let pointsTri = 
[
    // x      y       z
    [0.0,    0.75,   0.0],
    [-1.0,   0.0,    0.0],
    [1.0,   -1.0,   0.0]
];

const theta = 0.01;

const rotX = 
[
    [1, 0, 0],
    [0, Math.cos(theta), -Math.sin(theta)],
    [0, Math.sin(theta),  Math.cos(theta)]
];

const rotY = 
[
    [Math.cos(theta), 0, Math.sin(theta)],
    [0, 1, 0],
    [-Math.sin(theta), 0, Math.cos(theta)]
];

const rotZ = 
[
    [Math.cos(theta), -Math.sin(theta), 0],
    [Math.sin(theta), Math.cos(theta), 0],
    [0, 0, 1]
];


const color_1 = [255, 0, 0, 255];
const color_2 = [0, 255, 0, 255];
const color_3 = [0, 0, 255, 255];
const quality = 1000;

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
        
        render2DLine_noBlend(imageData, w, h, midpoint, padX, padY, pointsTri[0], pointsTri[1], color_1, quality);
        render2DLine_noBlend(imageData, w, h, midpoint, padX, padY, pointsTri[1], pointsTri[2], color_2, quality);
        render2DLine_noBlend(imageData, w, h, midpoint, padX, padY, pointsTri[2], pointsTri[0], color_3, quality);
        
        ctx.putImageData(imageData, 0, 0);
        lastTime = currentTime - (deltaTime % frameInterval);
    }
    
    requestAnimationFrame(render);
}

requestAnimationFrame(render);