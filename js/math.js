"use strict";
//plik zawierający przydatne funkcje matematyczne
/*jshint globalstrict: true*/
/*global window, THREE */

var trojkatPascala = [
  [1],
  [1, 1],
  [1, 2, 1],
  [1, 3, 3, 1]
];

function dwumian(n, k) {
  var s, wiersz;
  while (n >= trojkatPascala.length) {
    s = trojkatPascala.length;
    wiersz = [1];
    for (var i = 1; i < s; i++) {
      wiersz[i] = trojkatPascala[s - 1][i - 1] + trojkatPascala[s - 1][i];
    }
    wiersz[s] = 1;
    trojkatPascala.push(wiersz);
  }
  return trojkatPascala[n][k];
}

//wyznacz punkty posredni do wyznacznia wektorów A, B, C i D
//umożliwiających wyznacznie krzywizny
function wyznaczPunkty(inArray, u) {
  var temp = [inArray],
    temp2;
  for (var i = 0; i < inArray.length - 3; i++) {
    temp2 = [];
    for (var j = 0; j < temp[i].length - 1; j++) {
      temp2[j] = new THREE.Vector2();
      temp2[j].x = temp[i][j].x + (temp[i][j + 1].x - temp[i][j].x) * u;
      temp2[j].y = temp[i][j].y + (temp[i][j + 1].y - temp[i][j].y) * u;
    }
    temp.push(temp2);
  }
  return temp[temp.length - 1];
}

// Wyznaczenie krzywizny na podstawie algorytmu de Casteljau
function wyznaczKrzywizne(curve, ctrlPoint, count, deg) {
  var rGeometry = new THREE.Geometry();
  for (var i = 0; i <= count; i++) {
    var u = i / count;
    var Vx = [],
      Vy = [],
      temp = [];
    if (deg === 2) {
      for (var j = 0; j < ctrlPoint.length; j++) {
        Vx[j] = ctrlPoint[j].x;
        Vy[j] = ctrlPoint[j].y;
      }
    } else {
      temp = wyznaczPunkty(ctrlPoint, u);
      for (var k = 0; k < temp.length; k++) {
        Vx[k] = temp[k].x;
        Vy[k] = temp[k].y;
      }
    }

    var vecA = new THREE.Vector3(Vx[1] - Vx[0], Vy[1] - Vy[0]);
    var vecB = new THREE.Vector3(Vx[2] - Vx[1], Vy[2] - Vy[1]);

    Vx.push(Vx[0] + vecA.x * u);
    Vx.push(Vx[1] + vecB.x * u);
    Vy.push(Vy[0] + vecA.y * u);
    Vy.push(Vy[1] + vecB.y * u);

    var vecC = new THREE.Vector3(Vx[Vx.length - 1] - Vx[Vx.length - 2], Vy[Vy.length - 1] - Vy[Vy.length - 2]);
    var vecD = new THREE.Vector3().subVectors(vecB, vecA);

    var p1 = new THREE.Vector3(deg * vecC.x, deg * vecC.y);
    var p2 = new THREE.Vector3(deg * (deg - 1) * vecD.x, deg * (deg - 1) * vecD.y);

    var curv = (p1.x * p2.y - p1.y * p2.x) / Math.sqrt(Math.pow(p1.length(), 3));
    console.log("u: " + u + " k: " + curv);

    var temp = curve.getPoint(u);
    rGeometry.vertices.push(new THREE.Vector3(temp.x, temp.y));
    var temp2 = new THREE.Vector3();
    temp2.x = temp.x + p1.y * curv / p1.length();
    temp2.y = temp.y - p1.x * curv / p1.length();
    rGeometry.vertices.push(temp2);
  }
  return rGeometry;
}

function copyPoints(points) {
  var result = [];
  for (var i = 0; i < points.length; i++) {
    result.push(new THREE.Vector3(points[i].position.x, points[i].position.y, points[i].position.z));
  }
  return result;
}

//wyznacznie otoczki
function iloczyn(o, a, b) {
  return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

//algorytm Andrew's monotone chain
function wyznaczOtoczke(punkty) {
  punkty.sort(function (a, b) {
    return a.x == b.x ? a.y - b.y : a.x - b.x;
  });

  var gora = [];
  for (var i = 0; i < punkty.length; i++) {
    while (gora.length >= 2 && iloczyn(gora[gora.length - 2], gora[gora.length - 1], punkty[i]) >= 0) {
      gora.pop();
    }
    gora.push(punkty[i]);
  }

  var dol = [];
  for (var i = punkty.length - 1; i >= 0; i--) {
    while (dol.length >= 2 && iloczyn(dol[dol.length - 2], dol[dol.length - 1], punkty[i]) >= 0) {
      dol.pop();
    }
    dol.push(punkty[i]);
  }

  gora.pop();
  dol.pop();
  return gora.concat(dol);
}

function rysujOtoczke(points) {
  try {
    scene.remove(scene.getObjectByName('otoczka'));
  } catch (e) {
    console.log(e.message);
  }
  var otoczka = wyznaczOtoczke(points);
  var oGeometry = new THREE.Geometry();
  for (var i = 0; i < otoczka.length; i++) {
    otoczka[i].z = -10;
    oGeometry.vertices.push(otoczka[i]);
    if (i < otoczka.length - 2)
      oGeometry.faces.push(new THREE.Face3(0, i + 1, i + 2));
  }
  var oMaterial = new THREE.MeshBasicMaterial({
    color: "#c5ebfb",
    side: THREE.DoubleSide,
    transparent: true,
  });
  var oMesh = new THREE.Mesh(oGeometry, oMaterial);
  oMesh.name = "otoczka";
  scene.add(oMesh);
}

//algorytm Boehma do wyznaczania otoczki wypuklej
function boehmAlgoritm(controlPoints, nurbs, degree) {
  try {
    scene.remove(scene.getObjectByName('otoczka'));
    group = [];
  } catch (e) {
    console.log(e.message);
  }

  var geometry, mesh, arr;
  var processed = [],
    processed2 = [];
  var temp = new THREE.Vector3();
  var group = new THREE.Group();
  group.name = "otoczka";
  scene.add(group);

  var material = new THREE.MeshBasicMaterial({
    color: "#c5ebfb",
    side: THREE.DoubleSide,
    transparent: true,
  });
  if (degree == 2) {
    geometry = new THREE.Geometry();
    if (controlPoints.length == 3) {
      for (var i = 0; i < 3; i++) {
        controlPoints[i].position.z = -10;
        geometry.vertices.push(controlPoints[i].position);
      }
      geometry.faces.push(new THREE.Face3(0, 1, 2));
    } else if (controlPoints.length > 3) {
      //wyznaczenie wierzcholkow Beziera
      for (var i = 1; i <= controlPoints.length - 3; i++) {
        temp.subVectors(controlPoints[i + 1].position, controlPoints[i].position);
        temp.divideScalar(2);
        processed.push(new THREE.Vector3().addVectors(controlPoints[i].position, temp));
      }

      //definiowanie geometrii
      geometry.vertices.push(controlPoints[0].position, controlPoints[1].position);
      for (var i = 0; i < processed.length; i++) {
        geometry.vertices.push(processed[i], controlPoints[i + 2].position);
        geometry.faces.push(new THREE.Face3(2 * i, 2 * i + 1, 2 * i + 2));
      }
      geometry.vertices.push(controlPoints[controlPoints.length - 1].position);
      var l = geometry.vertices.length - 1;
      geometry.faces.push(new THREE.Face3(l - 2, l - 1, l));
      for (var i = 0; i <= l; i++) {
        geometry.vertices[i].z = -10;
      }
    }
    mesh = new THREE.Mesh(geometry, material);
    group.add(mesh);
  } else if (degree == 3) {
    geometry = new THREE.Geometry();
    if (controlPoints.length == 4) {
      for (var i = 0; i < controlPoints.length; i++) {
        controlPoints[i].position.z = -10;
        geometry.vertices.push(controlPoints[i].position);
      }
      var otoczka = wyznaczOtoczke(geometry.vertices);
      for (var i = 0; i < otoczka.length; i++) {
        geometry.vertices[i] = otoczka[i];
        if (i < otoczka.length - 2)
          geometry.faces.push(new THREE.Face3(0, i + 1, i + 2));
      }
      //rysowanie otoczki
      mesh = new THREE.Mesh(geometry, material);
      group.add(mesh);
    } else if (controlPoints.length > 4) {

      //wyznaczenie wierzcholkow Beziera dla krzywej 3 stopnia
      temp.subVectors(controlPoints[2].position, controlPoints[1].position);
      temp.divideScalar(2);
      processed.push(new THREE.Vector3().addVectors(controlPoints[1].position, temp));
      for (var i = 2; i < controlPoints.length - 3; i++) {
        temp.subVectors(controlPoints[i + 1].position, controlPoints[i].position);
        temp.divideScalar(3);
        processed.push(new THREE.Vector3().addVectors(controlPoints[i].position, temp));
        processed.push(new THREE.Vector3().subVectors(controlPoints[i + 1].position, temp));
      }
      temp.subVectors(controlPoints[controlPoints.length - 2].position, controlPoints[controlPoints.length - 3].position);
      temp.divideScalar(2);
      processed.push(new THREE.Vector3().addVectors(controlPoints[controlPoints.length - 3].position, temp));

      //wyznaczone boki dzielimy na pol
      for (var i = 0; i < processed.length; i += 2) {
        temp.subVectors(processed[i + 1], processed[i]);
        temp.divideScalar(2);
        processed2.push(new THREE.Vector3().addVectors(processed[i], temp));
      }
      //definiowanie geometrii
      geometry.vertices.push(controlPoints[0].position, controlPoints[1].position);
      for (var i = 0; i < processed2.length; i++) {
        geometry.vertices.push(processed[2 * i], processed2[i], processed[2 * i + 1]);
      }
      var l = controlPoints.length;
      geometry.vertices.push(controlPoints[l - 2].position, controlPoints[l - 1].position);
      l = geometry.vertices.length;
      for (var i = 0; i < l; i++) {
        geometry.vertices[i].z = -10;
        if (i % 3 === 0 && i !== 0) {
          arr = new THREE.Geometry();
          arr.vertices.push(geometry.vertices[i - 3], geometry.vertices[i - 2], geometry.vertices[i - 1], geometry.vertices[i]);
          var otoczka = wyznaczOtoczke(arr.vertices);
          for (var j = 0; j < otoczka.length; j++) {
            arr.vertices[j] = otoczka[j];
            if (j < otoczka.length - 2)
              arr.faces.push(new THREE.Face3(0, j + 1, j + 2));
          }
          //rysowanie otoczki
          mesh = new THREE.Mesh(arr, material);
          group.add(mesh);
        }
      }
    }
  }
}

function BezierSurface(controlPoints, n, m, weights) {
  THREE.Geometry.call(this);
  this.type = "BezierSurface";
  this.parameters = {
    uDegree: n,
    vDegree: m,
    controlPoints: controlPoints,
    weights: weights
  };

  this.controlPoints = controlPoints;
  this.weights = weights;
  var segments = 3;
  var verts = this.vertices;
  var faces = this.faces;

  createVertices();
  createFaces();

  function createVertices() {
    for (var v = 0; v / (m * segments) <= 1; v++) {
      for (var u = 0; u / (n * segments) <= 1; u++) {
        //tutaj wyznaczenie punktow dla u i v
        //czyli wywoalnie funkcji P(u,v)
        verts.push(bezier(u / (n * segments), v / (m * segments), n, m));
      }
    }
  }

  function createFaces() {
    var N = n * segments + 1,
      M = m * segments + 1;
    for (var i = 0; i < M - 1; i++) {
      for (var j = 0; j < N - 1; j++) {
        var a = i * N + j;
        var b = i * N + j + 1;
        var c = (i + 1) * N + j + 1;
        var d = (i + 1) * N + j;

        faces.push(new THREE.Face3(a, b, d));
        faces.push(new THREE.Face3(b, c, d));
      }
    }
  }

  function bezier(u, v, uDegree, vDegree) {
    var sumX = 0,
      sumY = 0,
      sumZ = 0,
      sumW = 0,
      bernU, bernV, cp, w;
    if (weights === undefined) {
      for (var i = 0; i <= vDegree; i++) {
        bernV = bernstein(vDegree, i, v);
        for (var j = 0; j <= uDegree; j++) {
          //wyznacz Bernsteina
          bernU = bernstein(uDegree, j, u);
          //liczba kolumn to uDegree+1
          cp = controlPoints[j + i * (uDegree + 1)].position;
          sumX += cp.x * bernU * bernV;
          sumY += cp.y * bernU * bernV;
          sumZ += cp.z * bernU * bernV;
        }
      }
      return new THREE.Vector3(sumX, sumY, sumZ);
    } else {
      for (var i = 0; i <= vDegree; i++) {
        bernV = bernstein(vDegree, i, v);
        for (var j = 0; j <= uDegree; j++) {
          //wyznacz Bernsteina
          bernU = bernstein(uDegree, j, u);
          //liczba kolumn to uDegree+1
          cp = controlPoints[j + i * (uDegree + 1)].position;
          w = weights[j + i * (uDegree + 1)];
          sumX += cp.x * bernU * bernV * w;
          sumY += cp.y * bernU * bernV * w;
          sumZ += cp.z * bernU * bernV * w;
          sumW += bernU * bernV * w;
        }
      }
      return new THREE.Vector3(sumX / sumW, sumY / sumW, sumZ / sumW);
    }
  }
}

function bernstein(N, i, u) {
  var k = 1 - u;
  return dwumian(N, i) * Math.pow(k, N - i) * Math.pow(u, i);
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
BezierCurve.prototype.getPoint = function (t) {
  var val = 0,
    sumaX = 0,
    sumaY = 0,
    suma = 0;
  var k = (1 - t);
  if (this.weight === null) {
    for (var i = 0; i <= this.n; i++) {
      val = dwumian(this.n, i) * Math.pow(k, (this.n - i)) * Math.pow(t, i);
      sumaX += val * this.geometry[i * 3];
      sumaY += val * this.geometry[i * 3 + 1];
    }
    return new THREE.Vector3(sumaX, sumaY, 0);
  } else {
    for (var i = 0; i <= this.n; i++) {
      val = dwumian(this.n, i) * Math.pow(k, (this.n - i)) * Math.pow(t, i);
      suma += val * this.weight[i];
      sumaX += val * this.geometry[i * 3] * this.weight[i];
      sumaY += val * this.geometry[i * 3 + 1] * this.weight[i];
    }
    return new THREE.Vector3(sumaX / suma, sumaY / suma, 0);
  }
};

function BernsteinCurve(stopien, numer) {
  this.n = stopien;
  this.k = numer;
}
BernsteinCurve.prototype = Object.create(THREE.Curve.prototype);
BernsteinCurve.prototype.constructor = BernsteinCurve;
BernsteinCurve.prototype.getPoint = function (t) {
  var k = (1 - t);
  var val = dwumian(this.n, this.k) * Math.pow(k, (this.n - this.k)) * Math.pow(t, this.k);
  return new THREE.Vector3(t * dim + (paddingH - dim), val * dim - paddingV, 0);
};
