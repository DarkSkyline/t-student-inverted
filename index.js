function betacf(x, a, b) {
    var fpmin = 1e-30;
    var m = 1;
    var qab = a + b;
    var qap = a + 1;
    var qam = a - 1;
    var c = 1;
    var d = 1 - qab * x / qap;
    var m2, aa, del, h;
  
    // These q's will be used in factors that occur in the coefficients
    if (Math.abs(d) < fpmin)
      d = fpmin;
    d = 1 / d;
    h = d;
  
    for (; m <= 100; m++) {
      m2 = 2 * m;
      aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      // One step (the even one) of the recurrence
      d = 1 + aa * d;
      if (Math.abs(d) < fpmin)
        d = fpmin;
      c = 1 + aa / c;
      if (Math.abs(c) < fpmin)
        c = fpmin;
      d = 1 / d;
      h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      // Next step of the recurrence (the odd one)
      d = 1 + aa * d;
      if (Math.abs(d) < fpmin)
        d = fpmin;
      c = 1 + aa / c;
      if (Math.abs(c) < fpmin)
        c = fpmin;
      d = 1 / d;
      del = d * c;
      h *= del;
      if (Math.abs(del - 1.0) < 3e-7)
        break;
    }
  
    return h;
  };

function gammaln(x) {
    var j = 0;
    var cof = [
      76.18009172947146, -86.50532032941677, 24.01409824083091,
      -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5
    ];
    var ser = 1.000000000190015;
    var xx, y, tmp;
    tmp = (y = xx = x) + 5.5;
    tmp -= (xx + 0.5) * Math.log(tmp);
    for (; j < 6; j++)
      ser += cof[j] / ++y;
    return Math.log(2.5066282746310005 * ser / xx) - tmp;
  };

  function ibeta(x, a, b) {
    // Factors in front of the continued fraction.
    var bt = (x === 0 || x === 1) ?  0 :
      Math.exp(gammaln(a + b) - gammaln(a) -
               gammaln(b) + a * Math.log(x) + b *
               Math.log(1 - x));
    if (x < 0 || x > 1)
      return false;
    if (x < (a + 1) / (a + b + 2))
      // Use continued fraction directly.
      return bt * betacf(x, a, b) / a;
    // else use continued fraction after making the symmetry transformation.
    return 1 - bt * betacf(1 - x, b, a) / b;
  };

// Returns the inverse of the incomplete beta function
function ibetainv(p, a, b) {
    var EPS = 1e-8;
    var a1 = a - 1;
    var b1 = b - 1;
    var j = 0;
    var lna, lnb, pp, t, u, err, x, al, h, w, afac;
    if (p <= 0)
      return 0;
    if (p >= 1)
      return 1;
    if (a >= 1 && b >= 1) {
      pp = (p < 0.5) ? p : 1 - p;
      t = Math.sqrt(-2 * Math.log(pp));
      x = (2.30753 + t * 0.27061) / (1 + t* (0.99229 + t * 0.04481)) - t;
      if (p < 0.5)
        x = -x;
      al = (x * x - 3) / 6;
      h = 2 / (1 / (2 * a - 1)  + 1 / (2 * b - 1));
      w = (x * Math.sqrt(al + h) / h) - (1 / (2 * b - 1) - 1 / (2 * a - 1)) *
          (al + 5 / 6 - 2 / (3 * h));
      x = a / (a + b * Math.exp(2 * w));
    } else {
      lna = Math.log(a / (a + b));
      lnb = Math.log(b / (a + b));
      t = Math.exp(a * lna) / a;
      u = Math.exp(b * lnb) / b;
      w = t + u;
      if (p < t / w)
        x = Math.pow(a * w * p, 1 / a);
      else
        x = 1 - Math.pow(b * w * (1 - p), 1 / b);
    }
    afac = -gammaln(a) - gammaln(b) + gammaln(a + b);
    for(; j < 10; j++) {
      if (x === 0 || x === 1)
        return x;
      err = ibeta(x, a, b) - p;
      t = Math.exp(a1 * Math.log(x) + b1 * Math.log(1 - x) + afac);
      u = err / t;
      x -= (t = u / (1 - 0.5 * Math.min(1, u * (a1 / x - b1 / (1 - x)))));
      if (x <= 0)
        x = 0.5 * (x + t);
      if (x >= 1)
        x = 0.5 * (x + t + 1);
      if (Math.abs(t) < EPS * x && j > 0)
        break;
    }
    return x;
  };
  
  function inv(p, dof, twoTailed = true, inverted = true) {
    p = twoTailed ? p / 2 : p;
    var x = ibetainv(2 * Math.min(p, 1 - p), 0.5 * dof, 0.5);
    x = Math.sqrt(dof * (1 - x) / x);
    var result = (p > 0.5) ? x : -x;
    result = inverted ? result * -1 : result;
    return result;
  };

for(var i = 1; i <= 100; i++){
    var valor = inv(0.04552, i)
    console.log(i + ' -> '+ valor.toFixed(2))    
}
