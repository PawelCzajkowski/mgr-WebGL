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
