///<reference path="three\three.d.ts"/>
///<reference path="three\three-orbitcontrols.d.ts"/>

/* solver code based off of: Jos Stam's work
    http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf
    3d extension based off of:
    https://github.com/BlainMaguire/3dfluid
*/

var stepping = false;
var applyRightForce = false;
var addFluid = false;
var timestep = false;
var fluidMode = true;
var vectorMode = false;
var hideFluid = false;
var hideArrows = false;
var maxParticlesPerBlock = 150;
var tex = THREE.ImageUtils.loadTexture(
    "cloud_1_512.png"
);

class FluidDynamicsSolver {
    dens: number[][][]; //x, y, z, from 0 to n + 2
    dens_prev: number[][][];
    v: number[][][];
    u: number[][][];
    w: number[][][];
    v_prev: number[][][];
    u_prev: number[][][];
    w_prev: number[][][];
    dt: number = 0.4; //fixed timestep
    diff: number = 0.0;
    visc: number = 0.0;
    lin_solver_times: number = 1;

    constructor(public n: number) {
        this.dens = [];
        this.dens_prev = [];
        this.v = [];
        this.u = [];
        this.w = [];
        this.v_prev = [];
        this.u_prev = [];
        this.w_prev = [];
        for (let i: number = 0; i < n + 2; i++) {
            this.dens[i] = [];
            this.dens_prev[i] = [];
            this.v[i] = [];
            this.u[i] = [];
            this.w[i] = [];
            this.v_prev[i] = [];
            this.u_prev[i] = [];
            this.w_prev[i] = [];
            for (let j: number = 0; j < n + 2; j++) {
                this.dens[i][j] = [];
                this.dens_prev[i][j] = [];
                this.v[i][j] = [];
                this.u[i][j] = [];
                this.w[i][j] = [];
                this.v_prev[i][j] = [];
                this.u_prev[i][j] = [];
                this.w_prev[i][j] = [];
                for (let k: number = 0; k < n + 2; k++) {
                    this.dens[i][j][k] = 0;
                    this.dens_prev[i][j][k] = 0;
                    this.v[i][j][k] = 0;
                    this.u[i][j][k] = 0;
                    this.w[i][j][k] = 0;
                    this.v_prev[i][j][k] = 0;
                    this.u_prev[i][j][k] = 0;
                    this.w_prev[i][j][k] = 0;
                }
            }
        }
    }

    add_source(N: number, x: number[][][], s: number[][][], dt: number) {
        let i, j, k;
        for (i = 1; i < N + 1; i++)
            for (j = 1; j < N + 1; j++)
                for (k = 1; k < N + 1; k++)
                    x[i][j][k] += dt * s[i][j][k];
    }

    set_bnd(N: number, b: number, x: number[][][]) {

        let i, j;
        for (i = 1; i <= N; i++) {
            for (j = 1; j <= N; j++) {
                x[0][i][j] = b == 1 ? -x[1][i][j] : x[1][i][j];
                x[N + 1][i][j] = b == 1 ? -x[N][i][j] : x[N][i][j];

                x[i][0][j] = b == 2 ? -x[i][1][j] : x[i][1][j];
                x[i][N + 1][j] = b == 2 ? -x[i][N][j] : x[i][N][j];

                x[i][j][0] = b == 3 ? -x[i][j][1] : x[i][j][1];
                x[i][j][N + 1] = b == 3 ? -x[i][j][N] : x[i][j][N];

            }
        }

        //corners
        x[0][0][0] = (1 / 3) * (x[1][0][0] + x[0][1][0] + x[0][0][1]);
        x[0][N + 1][0] = (1 / 3) * (x[1][N + 1][0] + x[0][N][0] + x[0][N + 1][1]);

        x[N + 1][0][0] = (1 / 3) * (x[N][0][0] + x[N + 1][1][0] + x[N + 1][0][1]);
        x[N + 1][N + 1][0] = (1 / 3) * (x[N][N + 1][0] + x[N + 1][N][0] + x[N + 1][N + 1][1]);

        x[0][0][N + 1] = (1 / 3) * (x[1][0][N + 1] + x[0][1][N + 1] + x[0][0][N]);
        x[0][N + 1][N + 1] = (1 / 3) * (x[1][N + 1][N + 1] + x[0][N][N + 1] + x[0][N + 1][N]);

        x[N + 1][0][N + 1] = (1 / 3) * (x[N][0][N + 1] + x[N + 1][1][N + 1] + x[N + 1][0][N]);
        x[N + 1][N + 1][N + 1] = (1 / 3) * (x[N][N + 1][N + 1] + x[N + 1][N][N + 1] + x[N + 1][N + 1][N]);

    }

    lin_solve(N: number, b: number, x: number[][][], x0: number[][][], a: number, c: number) {

        for (let t = 0; t < this.lin_solver_times; t++) {
            for (let i = 1; i <= N; i++) {
                for (let j = 1; j <= N; j++) {
                    for (let k = 1; k <= N; k++) {
                        x[i][j][k] = (x0[i][j][k] +
                            a * (x[i - 1][j][k] +
                                x[i + 1][j][k] +
                                x[i][j - 1][k] +
                                x[i][j + 1][k] +
                                x[i][j][k - 1] +
                                x[i][j][k + 1])
                        ) / c; 
                    }
                }
            }

            this.set_bnd(N, b, x);
        }
    }

    diffuse(N: number, b: number, x: number[][][], x0: number[][][], diff: number, dt: number) {
        let a: number = dt * diff * N * N * N;

        this.lin_solve(N, b, x, x0, a, 1 + 6 * a);
    }

    advect(N: number, b: number, d: number[][][], d0: number[][][], u: number[][][], v: number[][][], w: number[][][], dt: number) {
        let i, j, k, i0, j0, k0, i1, j1, k1; //int
        let x, y, z, s0, t0, u0, s1, t1, u1, dt0; //float

        dt0 = dt * N;
        for (i = 1; i <= N; i++) {
            for (j = 1; j <= N; j++) {
                for (k = 1; k <= N; k++) {
                    x = i - dt0 * u[i][j][k]; y = j - dt0 * v[i][j][k]; z = k - dt0 * w[i][j][k];

                    if (x < 0.5) x = 0.5; if (x > N + 0.5) x = N + 0.5; i0 = Math.floor(x); i1 = i0 + 1;
                    if (y < 0.5) y = 0.5; if (y > N + 0.5) y = N + 0.5; j0 = Math.floor(y); j1 = j0 + 1;
                    if (z < 0.5) z = 0.5; if (z > N + 0.5) z = N + 0.5; k0 = Math.floor(z); k1 = k0 + 1;

                    s1 = x - i0; s0 = 1 - s1;
                    t1 = y - j0; t0 = 1 - t1;
                    u1 = z - k0; u0 = 1 - u1;

                    d[i][j][k] = s0 * (
                        t0 * u0 * d0[i0][j0][k0] +
                        t1 * u0 * d0[i0][j1][k0] +
                        t0 * u1 * d0[i0][j0][k1] +
                        t1 * u1 * d0[i0][j1][k1]
                    ) + s1 * (
                        t0 * u0 * d0[i1][j0][k0] +
                        t1 * u0 * d0[i1][j1][k0] +
                        t0 * u1 * d0[i1][j0][k1] +
                        t1 * u1 * d0[i1][j1][k1]); //TODO: probably not correct for 3d, double check
                }
            }
        }

        this.set_bnd(N, b, d);
    }

    project(N: number, u: number[][][], v: number[][][], w: number[][][], p: number[][][], div: number[][][]) {

        let i, j, k;

        for (i = 1; i <= N; i++) {
            for (j = 1; j <= N; j++) {
                for (k = 1; k <= N; k++) {
                    div[i][j][k] = (-1.0 / 3.0) * ((u[i + 1][j][k] - u[i - 1][j][k]) / N + (v[i][j + 1][k] - v[i][j - 1][k]) / N + (w[i][j][k + 1] - w[i][j][k - 1]) / N);
                    p[i][j][k] = 0;
                }
            }
        }

        this.set_bnd(N, 0, div); this.set_bnd(N, 0, p);

        this.lin_solve(N, 0, p, div, 1, 6);

        for (i = 1; i <= N; i++) {
            for (j = 1; j <= N; j++) {
                for (k = 1; k <= N; k++) {
                    u[i][j][k] -= 0.5 * N * (p[i + 1][j][k] - p[i - 1][j][k]);
                    v[i][j][k] -= 0.5 * N * (p[i][j + 1][k] - p[i][j - 1][k]);
                    w[i][j][k] -= 0.5 * N * (p[i][j][k + 1] - p[i][j][k - 1]);
                }
            }
        }

        this.set_bnd(N, 1, u); this.set_bnd(N, 2, v); this.set_bnd(N, 3, w);

    }

    dens_step(N: number, x: number[][][], x0: number[][][], u: number[][][], v: number[][][], w: number[][][], diff: number, dt: number) {
        this.add_source(N, x, x0, dt);
        //SWAP( x0, x);
        this.diffuse(N, 0, x0, x, diff, dt); //NOTE, x0, x are swapped from Jos Stam's paper to simulate swaping pointers
        //SWAP( x0, x);
        this.advect(N, 0, x, x0, u, v, w, dt);
    }

    vel_step(N: number, u: number[][][], v: number[][][], w: number[][][], u0: number[][][], v0: number[][][], w0: number[][][], visc: number, dt: number) {
        this.add_source(N, u, u0, dt); this.add_source(N, v, v0, dt); this.add_source(N, w, w0, dt);
        //SWAP (u0, u);
        this.diffuse(N, 1, u0, u, visc, dt);
        //SWAP (v0, v);
        this.diffuse(N, 2, v0, v, visc, dt);
        //SWAP (w0, w);
        this.diffuse(N, 3, w0, w, visc, dt);
        this.project(N, u0, v0, w0, u, v); //TODO: example code doesn't use w0
        //SWAP (u0, u); SWAP (v0, v); SWAP (w0, w); //back to normal
        this.advect(N, 1, u, u0, u0, v0, w0, dt);
        this.advect(N, 2, v, v0, u0, v0, w0, dt);
        this.advect(N, 3, w, w0, u0, v0, w0, dt);
        this.project(N, u, v, w, u0, v0); //TODO: similar hariy problem
    }
}

class ParticleManager {
    particleCount: number;
    particles: THREE.Geometry;
    pMaterial: THREE.ParticleBasicMaterial;

    particleSystem: THREE.Points;

    constructor(public scene: THREE.Scene) {
        this.particleCount = 1800;
        this.particles = new THREE.Geometry();
        this.pMaterial = new THREE.PointsMaterial({
            color: 0xFFFFFF,
            size: 3,
            map: THREE.ImageUtils.loadTexture(
                "particle.png"
            ),
            blending: THREE.AdditiveBlending,
            transparent: true
        });

        for (var p = 0; p < this.particleCount; p++) {

            // create a particle with random
            // position values, -250 -> 250
            var pX = Math.random() * 500 - 250,
                pY = Math.random() * 500 - 250,
                pZ = Math.random() * 500 - 250,
                particle = new THREE.Vector3(
                    pX, pY, pZ
                );

            // add it to the geometry
            //this.particles.vertices.push(particle);
        }

        // create the particle system
        this.particleSystem = new THREE.Points(
            this.particles,
            this.pMaterial);

        // add it to the scene
        this.scene.add(this.particleSystem);
    }

    update() {
        this.particleSystem.rotation.x += 0.001;
    }

}

class ProjectWindow {
    scene: THREE.Scene;
    camera: THREE.Camera;
    renderer: THREE.Renderer;
    controls: THREE.OrbitControls;

    pManager: ParticleManager;
    cube: THREE.Mesh;
    grid: THREE.GridHelper;

    solver: FluidDynamicsSolver;

    cubeManager: THREE.Object3D;
    arrowManager: THREE.Object3D;

    fluidCubes: THREE.Mesh[];
    forceArrows: THREE.ArrowHelper[];

    constructor(element: HTMLElement) {
        this.scene = new THREE.Scene();
        this.camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);

        this.camera.translateX(300); this.camera.translateY(300); this.camera.translateZ(300);
        this.camera.lookAt(new THREE.Vector3(0, 0, 0));

        this.renderer = new THREE.WebGLRenderer();
        this.renderer.setSize(window.innerWidth, window.innerHeight);
        element.appendChild(this.renderer.domElement);

        this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableDamping = true;
        this.controls.dampingFactor = 0.25;
        this.controls.enableZoom = true;

        this.pManager = new ParticleManager(this.scene);

        var geometry = new THREE.BoxGeometry(1, 1, 1);
        var material = new THREE.MeshBasicMaterial({ color: 0x00ff00 });

        this.cube = new THREE.Mesh(geometry, material);

        this.grid = new THREE.GridHelper(200, 10, 0xff00ff, 0xff00ff);
        this.scene.add(this.grid);

        this.solver = new FluidDynamicsSolver(16);

        this.cubeManager = new THREE.Object3D();
        this.arrowManager = new THREE.Object3D();

        this.fluidCubes = [];
        this.forceArrows = [];
    }

    setupParticleSystem(n: number) {
        this.forceArrows = [];
        for (let i = 1; i <= n; i++) {
            for (let j = 1; j <= n; j++) {
                for (let k = 1; k <= n; k++) {
                    
                    let x = (i / n) * 500 - 250;
                    let y = (j / n) * 500 - 250;
                    let z = (k / n) * 500 - 250;

                    let c = new THREE.Color(0, 0, 0);
                    let geo = new THREE.BoxGeometry(500 / n, 500 / n, 500 / n);
                    let mat = new THREE.MeshBasicMaterial({
                        color: c.getHex(),
                        opacity: 0.05,
                        transparent: true
                    });

                    let cube = new THREE.Mesh(geo, mat);
                    cube.translateX(x); cube.translateY(y); cube.translateZ(z);
                    this.cubeManager.add(cube);
                    this.fluidCubes.push(cube);
                    this.scene.add(this.cubeManager);


                    //arrows

                    let from = new THREE.Vector3(x, y, z);
                    let to = new THREE.Vector3(x, y, z + 10);
                    var dir = to.clone().sub(from);
                    var len = dir.length();
                    var arrowHelper = new THREE.ArrowHelper(dir.normalize(), from, len, 0xff0000);
                    arrowHelper.visible = false;
                    this.arrowManager.add(arrowHelper);
                    this.forceArrows.push(arrowHelper);
                    this.scene.add(this.arrowManager);
                    
                    /*let geo = new THREE.Geometry();
                    let mat = new THREE.PointsMaterial({
                        color: 0xFFFFFF,
                        size: 10,
                        map: tex,
                        blending: THREE.AdditiveBlending,
                        transparent: true,
                        opacity: 0
                    });*/

                    //Was being used for initial particle placement
                    //let particleLoc = new THREE.Vector3(x, y, z);
                    //geo.vertices.push(particleLoc);

                   /* let part = new THREE.Points(geo, mat)
                    this.fluidParticleSystems.push(part);
                    this.scene.add(part);*/
                }
            }
        }

    }

    cleanPrev() {
        let n = this.solver.n;
        for (let i = 1; i <= n; i++) {
            for (let j = 1; j <= n; j++) {
                for (let k = 1; k <= n; k++) {
                    this.solver.dens_prev[i][j][k] = this.solver.u_prev[i][j][k] = this.solver.v_prev[i][j][k] = this.solver.w_prev[i][j][k] = 0;
                }
            }
        }
    }

    render() {

        var that = this;

        var n = that.solver.n;

        this.setupParticleSystem(n);
        this.cleanPrev();

        that.drawFluid();

        function step() {
                that.pManager.update();
                //if(stepping)
                //if(timestep)
                that.cleanPrev();

                if (hideArrows) {
                    that.arrowManager.traverse(function (object) {
                        object.visible = false;
                    });
                    that.cubeManager.traverse(function (object) {
                        object.visible = true;
                    });
                    hideArrows = false;
                }

                if (hideFluid) {
                    that.cubeManager.traverse(function (object) {
                        object.visible = false;
                    });
                    that.arrowManager.traverse(function (object) {
                        object.visible = true;
                    });
                    hideFluid = false;
                }



                if (addFluid)
                    that.solver.dens[n / 2][n / 2][n / 2] = 1500;

                if (applyRightForce)
                    that.solver.u[n / 2][2][n / 2] = 100;

                if (fluidMode)
                    that.drawFluid();

                if (vectorMode)
                    that.drawForces();


                that.renderer.render(that.scene, that.camera);

                window.requestAnimationFrame(step);
        }
        window.requestAnimationFrame(step);

    }

    drawFluid(): void {
        /*for (let cube of this.fluidCubes) {
            cube.material.dispose();
            cube.geometry.dispose();
            this.scene.remove(cube);
        }*/
        

        let n = this.solver.n;

        //this.solver.v_prev[n / 2][2][n / 2] = 1000;
        //this.solver.dens_prev[n/2][n/2][n/2] = 10000;

        this.solver.vel_step(this.solver.n, this.solver.u, this.solver.v, this.solver.w,
            this.solver.u_prev, this.solver.v_prev, this.solver.w_prev, this.solver.visc, this.solver.dt);
        this.solver.dens_step(this.solver.n, this.solver.dens, this.solver.dens_prev,
            this.solver.u, this.solver.v, this.solver.w, this.solver.diff, this.solver.dt);

       /* this.fluidCubes = [];*/

        let mat = new THREE.PointsMaterial({
            color: 0xFFFFFF,
            size: 10,
            map: tex,
            blending: THREE.NormalBlending,
            transparent: true,
            opacity: 0.4
        });

        let count = 0;
        for (let i = 1; i <= n; i++) {
            for (let j = 1; j <= n; j++) {
                for (let k = 1; k <= n; k++) {
                    let dens = Math.min(Math.max(this.solver.dens[i][j][k], 0), 1);

                    if (dens != 0) {
                        let x = (i / n) * 500 - 250;
                        let y = (j / n) * 500 - 250;
                        let z = (k / n) * 500 - 250;

                        let c = new THREE.Color(dens, dens, dens);
                        let cubeMat: THREE.MeshBasicMaterial = <THREE.MeshBasicMaterial>this.fluidCubes[count].material;
                        cubeMat.color.setHex(c.getHex());
                        cubeMat.needsUpdate = true;
                    }

                    //this.fluidParticleSystems[count].material.opacity = dens;
                    //this.fluidParticleSystems[count].material.needsUpdate = true;

                    //spawn random particles as a function of dens
                    //clear particle system


                   /* let geo = new THREE.Geometry();

                    let blockLength = 500 / n;

                    let x = (i / n) * 500 - 250;
                    let y = (j / n) * 500 - 250;
                    let z = (k / n) * 500 - 250;

                    let totalParticles = Math.floor(dens * maxParticlesPerBlock);
                    for (let l = 0; l < totalParticles; l++) {

                        let px = Math.random() * blockLength - (blockLength / 2) + x,
                            py = Math.random() * blockLength - (blockLength / 2) + y,
                            pz = Math.random() * blockLength - (blockLength / 2) + z;

                        let new_part = new THREE.Vector3(px, py, pz);
                        geo.vertices.push(new_part);
                    }
                    this.fluidParticleSystems[count].material.dispose();
                    this.fluidParticleSystems[count].geometry.dispose();
                    this.scene.remove(this.fluidParticleSystems[count]);

                    this.fluidParticleSystems[count] = new THREE.Points(geo, mat);
                    this.scene.add(this.fluidParticleSystems[count]);*/

                    count++;
                }
            }
        }
    }

    drawForces(): void {
        /*for (let arrow of this.forceArrows) {
            this.scene.remove(arrow);
        }*/

        let i, j, k;
        let x, y, z, h;
        let N = this.solver.n;
        let u = this.solver.u;
        let v = this.solver.v;
        let w = this.solver.w;


        h = 1 / N;

        this.solver.vel_step(this.solver.n, this.solver.u, this.solver.v, this.solver.w,
            this.solver.u_prev, this.solver.v_prev, this.solver.w_prev, this.solver.visc, this.solver.dt);
        this.solver.dens_step(this.solver.n, this.solver.dens, this.solver.dens_prev,
            this.solver.u, this.solver.v, this.solver.w, this.solver.diff, this.solver.dt);


        let n = this.solver.n;
        let count = 0;
        for (i = 1; i <= N; i++) {
            for (j = 1; j <= N; j++) {
                for (k = 1; k <= N; k++) {

                    let x = (i / n) * 500 - 250;
                    let y = (j / n) * 500 - 250;
                    let z = (k / n) * 500 - 250;
                    let from = new THREE.Vector3(x, y, z);
                    let to = new THREE.Vector3(x + u[i][j][k], y + v[i][j][k], z + w[i][j][k]);
                    var dir = to.clone().sub(from).normalize();
                    var len = dir.length() * 10;
                    let arrow = this.forceArrows[count];

                    arrow.position = from;
                    arrow.setDirection(dir);
                    arrow.setLength(len);
                    count++
                }
            }
        }
        


    }

    testSolverRunTime(): number {
        let start = new Date();
        this.solver.vel_step(this.solver.n, this.solver.u, this.solver.v, this.solver.w,
            this.solver.u_prev, this.solver.v_prev, this.solver.w_prev, this.solver.visc, this.solver.dt);
        this.solver.dens_step(this.solver.n, this.solver.dens, this.solver.dens_prev,
            this.solver.u, this.solver.v, this.solver.w, this.solver.diff, this.solver.dt);
        let end = new Date();
        return end.getTime() - start.getTime();
    }

}

window.onload = () => {
    let el = document.getElementById('content');

    let win = new ProjectWindow(el);
    console.log(win.testSolverRunTime());
    win.render();

};

window.onkeydown = function (e) {
    var key = e.keyCode ? e.keyCode : e.which;

    if (key == 38) //up key
        stepping = true;
    else if (key == 80) //p key
        applyRightForce = true;
    else if (key == 79) // o
        addFluid = true;
    else if (key == 76)
        timestep = true;
    else if (key == 86) // v
    {
        vectorMode = true;
        fluidMode = false;
        hideFluid = true;
    }
    else if (key == 70)// f
    {
        fluidMode = true
        vectorMode = false;
        hideArrows = true;
    }
}

window.onkeyup = function (e) {
    var key = e.keyCode ? e.keyCode : e.which;

    if (key == 38) //up key
        stepping = false;
    else if (key == 80) //right key
        applyRightForce = false;
    else if (key == 79) //o key
        addFluid = false;
    else if (key == 76) // l
        timestep = false;

}