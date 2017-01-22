/* @author Paweł Czajkowski
 * Copyrights © 2016 All rights reserved.
 */

"use strict";
//plik zawierający przydatne funkcje matematyczne

var trojkatPascala = [
  [1],
  [1, 1],
  [1, 2, 1],
  [1, 3, 3, 1]
];

var material = new THREE.MeshBasicMaterial({
  color: "#c5ebfb",
  side: THREE.DoubleSide,
  transparent: true,
});

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

function length(v) {
  return Math.sqrt(v.x * v.x);
}

function checkKnots(arr, j) {
  for (var i = 0; i < arr.length; i++) {
    if (i < j && arr[i] > arr[j]) arr[i] = arr[j];
    if (i > j && arr[i] < arr[j]) arr[i] = arr[j];
  }
}

function sklejPunkty(arr, i) {
  var t;
  if (arr[i] <= 0.05) {
    arr[i] = 0;
  } else if (arr[i] >= 0.95) {
    arr[i] = 1;
  } else if ((arr[i] - arr[i - 1] <= 0.05)) {
    arr[i] = arr[i - 1];
  } else if ((arr[i + 1] - arr[i]) <= 0.05) {
    arr[i] = arr[i + 1];
  }
  temp = arr[i] - arr[i - 1];
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

//wymierny algorytm de Casteljau
function wymiernyDeCasteljau(punkty, u, wagi) {
  var W = [],
    tempW = [];
  var P = [],
    tempP = [];
  var w, s, p, temp1, temp2;

  for (var i = 0; i < punkty.length; i++) {
    tempP.push(punkty[i].clone());
    if (wagi === undefined) {
      tempW.push(1);
    } else {
      tempW.push(wagi[i]);
    }
  }
  P.push(tempP);
  W.push(tempW);

  //wyznaczenie punktow posrednich
  var n = punkty.length - 1;
  for (var j = 1; j <= n; j++) {
    tempP = [];
    tempW = [];
    for (var i = 0; i <= n - j; i++) {
      w = (1 - u) * W[j - 1][i] + u * W[j - 1][i + 1];
      s = u * W[j - 1][i + 1] / w;
      temp1 = P[j - 1][i + 1].clone();
      temp1.multiplyScalar(s);
      s = (1 - u) * W[j - 1][i] / w;
      temp2 = P[j - 1][i].clone();
      temp2.multiplyScalar(s);
      p = new THREE.Vector3().addVectors(temp1, temp2);
      tempP.push(p);
      tempW.push(w);
    }
    P.push(tempP);
    W.push(tempW);
  }
  return P;
}

// Wyznaczenie krzywizny na podstawie algorytmu de Casteljau
function wyznaczKrzywizne(curve, ctrlPoint, count, deg, factor, offset, number) {
  var rGeometry = new THREE.Geometry();
  if (offset == undefined) offset = 0;
  if (number == undefined) number = 1;
  if (factor == undefined) factor = 1;
  for (var i = 0; i <= count; i++) {
    var u = i / count;

    //dla krzywych wymiernych
    if (geometry.type === "Geometry" || geometry.getAttribute("weight") === undefined) {
      var V = wymiernyDeCasteljau(ctrlPoint, u);
    } else {
      var V = wymiernyDeCasteljau(ctrlPoint, u, geometry.attributes.weight.array);
    }
    var vecA = new THREE.Vector3().subVectors(V[deg - 2][1], V[deg - 2][0]);
    var vecB = new THREE.Vector3().subVectors(V[deg - 2][2], V[deg - 2][1]);
    var vecC = new THREE.Vector3().subVectors(V[deg - 1][1], V[deg - 1][0]);
    var vecD = new THREE.Vector3().subVectors(vecB, vecA);

    var p1 = new THREE.Vector3(deg * vecC.x, deg * vecC.y);
    var p2 = new THREE.Vector3(deg * (deg - 1) * vecD.x, deg * (deg - 1) * vecD.y);

    var curv = (p1.x * p2.y - p1.y * p2.x) / Math.sqrt(Math.pow(p1.length(), 3));

    var temp = curve.getPoint(i / (count * number) + offset / number);
    rGeometry.vertices.push(new THREE.Vector3(temp.x, temp.y));
    var temp2 = new THREE.Vector3();

    temp2.x = temp.x + p1.y * curv * factor / p1.length();
    temp2.y = temp.y - p1.x * curv * factor / p1.length();

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
  var oMesh = new THREE.Mesh(oGeometry, material);
  oMesh.name = "otoczka";
  scene.add(oMesh);
}

//funkcja wykorzystywana w aplecie 8.39
function wyswietlOtoczke(controlPoints, degree) {
  try {
    scene.remove(scene.getObjectByName('otoczka'));
    group = [];
  } catch (e) {
    console.log(e.message);
  }
  var mesh;
  var group = new THREE.Group();
  group.name = "otoczka";
  scene.add(group);
  for (var i = 0; i < controlPoints.length - degree; i++) {
    otoczkaGeometry = new THREE.Geometry();
    for (var j = 0; j < degree + 1; j++) {
      controlPoints[i + j].z = -10;
      otoczkaGeometry.vertices.push(controlPoints[i + j]);
    }
    var otoczka = wyznaczOtoczke(otoczkaGeometry.vertices);
    for (var k = 0; k < otoczka.length; k++) {
      otoczkaGeometry.vertices[k] = otoczka[k];
      if (k < otoczka.length - 2)
        otoczkaGeometry.faces.push(new THREE.Face3(0, k + 1, k + 2));
    }
    //rysowanie otoczki
    mesh = new THREE.Mesh(otoczkaGeometry, material);
    group.add(mesh);
  }
}

//algorytm Boehma do wyznaczania otoczki wypuklej
function boehmAlgoritm(nurbs, draw) {
  try {
    scene.remove(scene.getObjectByName("otoczka"));
    group = [];
  } catch (e) {
    console.log(e.message);
  }

  var pGeometry = new THREE.CircleBufferGeometry(5, 16);
  var pMaterial = new THREE.MeshBasicMaterial({
    color: "black"
  });
  if (draw === undefined) draw = true;

  var geometry, mesh, arr;
  var processed = [],
    processed2 = [];
  var d = [];
  var A, B, C, D;
  var temp = new THREE.Vector3();

  var group = new THREE.Group();
  group.name = "otoczka";
  scene.add(group);

  if (nurbs.knots !== undefined || nurbs.knots.length > 0) {
    for (var i = 0; i < nurbs.knots.length - 1; i++) {
      d.push(nurbs.knots[i + 1] - nurbs.knots[i]);
    }
  }
  var material = new THREE.MeshBasicMaterial({
    color: "#c5ebfb",
    side: THREE.DoubleSide,
    transparent: true,
  });
  if (nurbs.degree == 2) {
    if (nurbs.controlPoints.length == 3) {
      geometry = new THREE.Geometry();
      for (var i = 0; i < 3; i++) {
        nurbs.controlPoints[i].z = -10;
        geometry.vertices.push(nurbs.controlPoints[i]);
      }
      for (var j = 0; j < 3; j++) {
        geometry.vertices[j].z = -10;
      }
      geometry.faces.push(new THREE.Face3(0, 1, 2));
      mesh = new THREE.Mesh(geometry, material);
      group.add(mesh);

    } else if (nurbs.controlPoints.length > 3) {
      //wyznaczenie wierzcholkow Beziera
      for (var i = 1; i <= nurbs.controlPoints.length - 3; i++) {
        A = nurbs.controlPoints[i + 1], B = nurbs.controlPoints[i];
        C = new THREE.Vector4();
        C.w = (A.w * d[i + 2] + B.w * d[i + 1]) / (d[i + 1] + d[i + 2]);
        C.x = (A.x * A.w * d[i + 2] + B.x * B.w * d[i + 1]) / (C.w * (d[i + 1] + d[i + 2]));
        C.y = (A.y * A.w * d[i + 2] + B.y * B.w * d[i + 1]) / (C.w * (d[i + 1] + d[i + 2]));
        processed.push(C);
      }

      //definiowanie geometrii
      geometry = new THREE.Geometry();
      geometry.vertices.push(nurbs.controlPoints[0], nurbs.controlPoints[1], processed[0]);
      for (var j = 0; j < 3; j++) {
        geometry.vertices[j].z = -10;
      }
      geometry.faces.push(new THREE.Face3(0, 1, 2));
      mesh = new THREE.Mesh(geometry, material);
      group.add(mesh);
      for (var i = 0; i < processed.length - 1; i++) {
        geometry = new THREE.Geometry();
        geometry.vertices.push(processed[i], nurbs.controlPoints[i + 2], processed[i + 1]);
        for (var j = 0; j < 3; j++) {
          geometry.vertices[j].z = -10;
        }
        geometry.faces.push(new THREE.Face3(0, 1, 2));
        mesh = new THREE.Mesh(geometry, material);

        var point = new THREE.Mesh(pGeometry, pMaterial);
        point.position.copy(processed[i]);
        point.position.z = 2;
        group.add(mesh);
        group.add(point);
      }
      geometry = new THREE.Geometry();
      geometry.vertices.push(processed[processed.length - 1]);
      geometry.vertices.push(nurbs.controlPoints[nurbs.controlPoints.length - 2]);
      geometry.vertices.push(nurbs.controlPoints[nurbs.controlPoints.length - 1]);
      for (var j = 0; j < 3; j++) {
        geometry.vertices[j].z = -10;
      }
      geometry.faces.push(new THREE.Face3(0, 1, 2));
      mesh = new THREE.Mesh(geometry, material);

      var point = new THREE.Mesh(pGeometry, pMaterial);
      point.position.copy(processed[processed.length - 1]);
      point.position.z = 2;
      group.add(mesh);
      group.add(point);
    }
  } else if (nurbs.degree == 3) {
    geometry = new THREE.Geometry();
    if (nurbs.controlPoints.length == 4) {
      for (var i = 0; i < nurbs.controlPoints.length; i++) {
        nurbs.controlPoints[i].z = -10;
        geometry.vertices.push(nurbs.controlPoints[i]);
      }
      if (draw) {
        var otoczka = wyznaczOtoczke(geometry.vertices);
        for (var i = 0; i < otoczka.length; i++) {
          geometry.vertices[i] = otoczka[i];
          if (i < otoczka.length - 2)
            geometry.faces.push(new THREE.Face3(0, i + 1, i + 2));
        }
      }

      //rysowanie otoczki
      mesh = new THREE.Mesh(geometry, material);
      group.add(mesh);
    } else if (nurbs.controlPoints.length > 4) {
      for (var i = 1; i < nurbs.controlPoints.length - 2; i++) {
        A = nurbs.controlPoints[i], B = nurbs.controlPoints[i + 1];
        C = new THREE.Vector4(), D = new THREE.Vector4();
        C.w = (A.w * (d[i + 2] + d[i + 3]) + B.w * d[i + 1]) / (d[i + 1] + d[i + 2] + d[i + 3]);
        D.w = (A.w * d[i + 3] + B.w * (d[i + 1] + d[i + 2])) / (d[i + 1] + d[i + 2] + d[i + 3]);
        C.x = (A.x * A.w * (d[i + 2] + d[i + 3]) + B.x * B.w * d[i + 1]) / (C.w * (d[i + 1] + d[i + 2] + d[i + 3]));
        C.y = (A.y * A.w * (d[i + 2] + d[i + 3]) + B.y * B.w * d[i + 1]) / (C.w * (d[i + 1] + d[i + 2] + d[i + 3]));
        D.x = (A.x * A.w * d[i + 3] + B.x * B.w * (d[i + 1] + d[i + 2])) / (D.w * (d[i + 1] + d[i + 2] + d[i + 3]));
        D.y = (A.y * A.w * d[i + 3] + B.y * B.w * (d[i + 1] + d[i + 2])) / (D.w * (d[i + 1] + d[i + 2] + d[i + 3]));
        processed.push(C);
        processed.push(D);

        var point = new THREE.Mesh(pGeometry, pMaterial);
        point.position.copy(C);
        group.add(point);

        point = new THREE.Mesh(pGeometry, pMaterial);
        point.position.copy(D);
        group.add(point);
      }

      //wyznaczone boki dzielimy na pol
      for (var i = 0; i < processed.length / 2 - 1; i++) {
        A = processed[2 * i + 1], B = processed[2 * (i + 1)];
        C = new THREE.Vector4();
        C.w = (A.w * (d[i + 4]) + B.w * d[i + 3]) / (d[i + 3] + d[i + 4]);
        C.x = (A.x * A.w * (d[i + 4]) + B.x * B.w * d[i + 3]) / (C.w * (d[i + 3] + d[i + 4]));
        C.y = (A.y * A.w * (d[i + 4]) + B.y * B.w * d[i + 3]) / (C.w * (d[i + 3] + d[i + 4]));
        processed2.push(C);

        var point = new THREE.Mesh(pGeometry, pMaterial);
        point.position.copy(C);
        point.position.z = 2;
        group.add(point);
      }

      //definiowanie geometrii
      geometry.vertices.push(nurbs.controlPoints[0], nurbs.controlPoints[1]);
      for (var i = 0; i < processed2.length; i++) {
        geometry.vertices.push(processed[2 * i + 1], processed2[i], processed[2 * i + 2]);
      }
      var l = nurbs.controlPoints.length;
      geometry.vertices.push(nurbs.controlPoints[l - 2], nurbs.controlPoints[l - 1]);
      l = geometry.vertices.length;
      for (var i = 0; i < l; i++) {
        geometry.vertices[i].z = -10;
        if (i % 3 === 0 && i !== 0) {
          arr = new THREE.Geometry();
          arr.vertices.push(geometry.vertices[i - 3], geometry.vertices[i - 2], geometry.vertices[i - 1], geometry.vertices[i]);
          if (draw) {
            var otoczka = wyznaczOtoczke(arr.vertices);
            for (var j = 0; j < otoczka.length; j++) {
              arr.vertices[j] = otoczka[j];
              if (j < otoczka.length - 2)
                arr.faces.push(new THREE.Face3(0, j + 1, j + 2));
            }
          }
          //rysowanie otoczki
          mesh = new THREE.Mesh(arr, material);
          group.add(mesh);

        }
      }
      //line segments
      var sGeometry = new THREE.Geometry();
      for (var i = 1; i < processed.length - 1; i++) {
        sGeometry.vertices.push(processed[i]);
      }
      var sMaterial = new THREE.LineBasicMaterial({
        color: "black",
        linewidth: 2
      });
      group.add(new THREE.LineSegments(sGeometry, sMaterial));
    }
  }
  if (!draw) {
    group.visible = false;
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
  var uSegments = Number(document.getElementById("uSegments").value);
  var vSegments = Number(document.getElementById("vSegments").value);
  var verts = this.vertices;
  var faces = this.faces;

  createVertices();
  createFaces();

  function createVertices() {
    for (var v = 0; v / (m * vSegments) <= 1; v++) {
      for (var u = 0; u / (n * uSegments) <= 1; u++) {
        //tutaj wyznaczenie punktow dla u i v
        //czyli wywoalnie funkcji P(u,v)
        verts.push(bezier(u / (n * uSegments), v / (m * vSegments), n, m));
      }
    }
  }

  function createFaces() {
    var N = n * uSegments + 1,
      M = m * vSegments + 1;
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
      val = dwumian(this.n, i) * Math.pow(k, (this.n - i)) * Math.pow(t, i) * this.weight[i];
      suma += val;
      sumaX += val * this.geometry[i * 3];
      sumaY += val * this.geometry[i * 3 + 1];
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
  return new THREE.Vector3(A.x + t * V.x, A.y + val * Vy.y, 0);
};
