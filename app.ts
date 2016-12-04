///<reference path="three\three.d.ts"/>
///<reference path="three\three-orbitcontrols.d.ts"/>


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
                particle = new THREE.Vertex(
                    pX, pY, pZ
                );

            // add it to the geometry
            this.particles.vertices.push(particle);
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

    constructor(element: HTMLElement) {
        this.scene = new THREE.Scene();
        this.camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
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
    }

    addCube() {
        //this.scene.add(this.cube);

        this.camera.position.z = 5;
    }

    render() {

        var that = this;
        setInterval(function () {
            that.pManager.update();
            that.renderer.render(that.scene, that.camera);
        }, 1000 / 60);


    }

}

window.onload = () => {
    let el = document.getElementById('content');

    let win = new ProjectWindow(el);
    win.addCube();
    win.render();

};