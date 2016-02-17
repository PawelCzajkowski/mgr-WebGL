//plik zawierający przydatne funkcje matematyczne
/*jshint globalstrict: true*/
"use strict";

var trojkatPascala = [
  [1],
  [1, 1],
  [1, 2, 1],
  [1, 3, 3, 1]
];

function dwumian(n, k) {
  while (n >= trojkatPascala.length) {
    var s = trojkatPascala.length;
    var wiersz = [1];
    for (var i = 1; i < s; i++ ) {
      wiersz[i] = trojkatPascala[s-1][i-1] + trojkatPascala[s-1][i];
    }
    wiersz[s]=1;
    trojkatPascala.push(wiersz);
  }
  return trojkatPascala[n][k];
}

function BezierSurface(controlPoints, n, m, weights) {
  THREE.Geometry.call(this);
  this.type = "BezierSurface";
  this.parameters = {
    uDegree : n,
    vDegree : m,
    controlPoints : controlPoints,
    weights : weights
  };

  this.controlPoints = controlPoints;
  this.weights = weights || 1;
  var segments = 3;
  var verts = this.vertices;
  var faces = this.faces;

  createVertices();
  createFaces();

  function createVertices() {
    for (var v=0; v/(m*segments) <= 1; v++){
      for (var u=0; u/(n*segments) <= 1; u++) {
        //tutaj wyznaczenie punktow dla u i v
        //czyli wywoalnie funkcji P(u,v)
        verts.push(bezier(u/(n*segments), v/(m*segments), n, m));
      }
    }
  }

  function createFaces() {
    var N = n*segments+1, M = m*segments+1;
    for (var i=0; i<M-1; i++) {
      for (var j=0; j<N-1; j++) {
        var a = i*N+j;
        var b = i*N+j+1;
        var c = (i+1)*N+j+1;
        var d = (i+1)*N+j;

        faces.push(new THREE.Face3(a,b,d));
        faces.push(new THREE.Face3(b,c,d));
      }
    }
  }

  function bezier(u, v, uDegree, vDegree) {
    var sumX = 0, sumY = 0, sumZ = 0;
    for (var i=0; i <= vDegree; i++) {
      var bernV = bernstein(vDegree, i, v);
      for (var j=0; j <= uDegree; j++) {
        //wyznacz Bernsteina
        var bernU = bernstein(uDegree, j, u);
        //liczba kolumn to uDegree+1
        var cp = controlPoints[j+i*(uDegree+1)].position;
        sumX += cp.x*bernU*bernV;
        sumY += cp.y*bernU*bernV;
        sumZ += cp.z*bernU*bernV;
      }
    }
    return new THREE.Vector3(sumX, sumY, sumZ);
  }

  function bernstein(N, i, u) {
    var k = 1-u;
    return dwumian(N, i)*Math.pow(k, N-i)*Math.pow(u, i);
  }
}

BezierSurface.prototype = Object.create(THREE.Geometry.prototype);
BezierSurface.prototype.constructor = BezierSurface;

function BezierCurve(stopien, tablica, waga) {
  this.n = stopien;
  this.geometry = tablica.position.array;
  if (typeof waga === "undefined") {
    this.weight = null;
  } else {
    this.weight = tablica.weight.array;
  }
}

BezierCurve.prototype = Object.create(THREE.Curve.prototype);
BezierCurve.prototype.constructor = BezierCurve;
BezierCurve.prototype.getPoint  = function(t) {
  var val = 0, sumaX = 0, sumaY = 0, suma = 0;
  var k = (1-t);
  if (this.weight === null)
  {
    for (var i=0; i<=this.n; i++) {
      val = dwumian(this.n,i)*Math.pow(k,(this.n-i))*Math.pow(t,i);
      sumaX += val*this.geometry[i*3];
      sumaY += val*this.geometry[i*3+1];
    }
    return new THREE.Vector3(sumaX, sumaY, 0);
  }
  else
  {
    for (var i=0; i<=this.n; i++) {
      val = dwumian(this.n,i)*Math.pow(k,(this.n-i))*Math.pow(t,i);
      suma += val*this.weight[i];
      sumaX += val*this.geometry[i*3]*this.weight[i];
      sumaY += val*this.geometry[i*3+1]*this.weight[i];
    }
    return new THREE.Vector3(sumaX/suma, sumaY/suma, 0);
  }
};

function BernsteinCurve(stopien, numer) {
  this.n = stopien;
  this.k = numer;
}
BernsteinCurve.prototype = Object.create(THREE.Curve.prototype);
BernsteinCurve.prototype.constructor = BernsteinCurve;
BernsteinCurve.prototype.getPoint = function(t) {
  var k = (1-t);
  var val = dwumian(this.n,this.k)*Math.pow(k,(this.n-this.k))*Math.pow(t,this.k);
  return new THREE.Vector3(t*dim+(paddingH-dim), val*dim-paddingV, 0);
};
