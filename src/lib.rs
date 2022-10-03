// Translated from the Cephes library 2.1 C source code, which has the
// following notice:
//   Cephes Math Library Release 2.1:  January, 1989
//   Copyright 1984, 1987, 1989 by Stephen L. Moshier

type R = f64;
const PI : R = std::f64::consts::PI;
const FRAC_PI_2 : R = std::f64::consts::FRAC_PI_2;

fn polevl(x:R,coef:&[R],n:usize)->R {
    let mut ans = coef[0];
    for &c in &coef[1..=n] {
	ans = ans*x + c;
    }
    ans
}

fn p1evl(x:R,coef:&[R],n:usize)->R {
    let mut ans = x + coef[0];
    for &c in &coef[1..n] {
	ans = ans*x + c;
    }
    ans
}

/// Evaluates the Fresnel integrals C(x)=∫cos(π/2 t²)dt and
/// S(x)=∫sin(π/2 t²)dt for t∈\[0,x\].  Returns (S(x),C(x)).
/// Based on the Cephes library function of the same name.
pub fn fresnl(x:R)->(R,R) {
    let xxa = x;
    
    let mut cc;
    let mut ss;

    // S(x) for small x
    let s_n : [R;6] = [
	-2.99181919401019853726E3,
	7.08840045257738576863E5,
	-6.29741486205862506537E7,
	2.54890880573376359104E9,
	-4.42979518059697779103E10,
	3.18016297876567817986E11,
    ].into();

    let s_d : [R;6] = [
	// 1.00000000000000000000E0,
	2.81376268889994315696E2,
	4.55847810806532581675E4,
	5.17343888770096400730E6,
	4.19320245898111231129E8,
	2.24411795645340920940E10,
	6.07366389490084639049E11,
    ];

    // C(x) for small x
    let c_n : [R;6] = [
	-4.98843114573573548651E-8,
	9.50428062829859605134E-6,
	-6.45191435683965050962E-4,
	1.88843319396703850064E-2,
	-2.05525900955013891793E-1,
	9.99999999999999998822E-1,
    ];

    let c_d : [R;7] = [
	3.99982968972495980367E-12,
	9.15439215774657478799E-10,
	1.25001862479598821474E-7,
	1.22262789024179030997E-5,
	8.68029542941784300606E-4,
	4.12142090722199792936E-2,
	1.00000000000000000118E0,
    ];

    // Auxiliary function f(x)
    let f_n : [R;10] = [
	4.21543555043677546506E-1,
	1.43407919780758885261E-1,
	1.15220955073585758835E-2,
	3.45017939782574027900E-4,
	4.63613749287867322088E-6,
	3.05568983790257605827E-8,
	1.02304514164907233465E-10,
	1.72010743268161828879E-13,
	1.34283276233062758925E-16,
	3.76329711269987889006E-20,
    ];

    let f_d : [R;10] = [
	/*  1.00000000000000000000E0, */
	7.51586398353378947175E-1,
	1.16888925859191382142E-1,
	6.44051526508858611005E-3,
	1.55934409164153020873E-4,
	1.84627567348930545870E-6,
	1.12699224763999035261E-8,
	3.60140029589371370404E-11,
	5.88754533621578410010E-14,
	4.52001434074129701496E-17,
	1.25443237090011264384E-20,
    ];

    /* Auxiliary function g(x) */
    let g_n : [R;11] = [
	5.04442073643383265887E-1,
	1.97102833525523411709E-1,
	1.87648584092575249293E-2,
	6.84079380915393090172E-4,
	1.15138826111884280931E-5,
	9.82852443688422223854E-8,
	4.45344415861750144738E-10,
	1.08268041139020870318E-12,
	1.37555460633261799868E-15,
	8.36354435630677421531E-19,
	1.86958710162783235106E-22,
    ];

    let g_d : [R;11] = [
	/*  1.00000000000000000000E0, */
	1.47495759925128324529E0,
	3.37748989120019970451E-1,
	2.53603741420338795122E-2,
	8.14679107184306179049E-4,
	1.27545075667729118702E-5,
	1.04314589657571990585E-7,
	4.60680728146520428211E-10,
	1.10273215066240270757E-12,
	1.38796531259578871258E-15,
	8.39158816283118707363E-19,
	1.86958710162783236342E-22,
    ];
    
    loop {
	if xxa.is_infinite() {
	    cc = 0.5;
	    ss = 0.5;
	    break;
	}

	let x = xxa.abs();
	let x2 = x * x;
	if x2 < 2.5625 {
	    let t = x2 * x2;
	    ss = x * x2 * polevl(t,&s_n,5) / p1evl(t,&s_d,6);
	    cc = x * polevl(t,&c_n,5) / polevl(t,&c_d,6);
	    break;
	}

	if x > 36974.0 {
	    // http://functions.wolfram.com/GammaBetaErf/FresnelC/06/02/
	    // http://functions.wolfram.com/GammaBetaErf/FresnelS/06/02/
	    cc = 0.5 + 1.0/(PI*x) * (PI*x*x/2.0).sin();
	    ss = 0.5 - 1.0/(PI*x) * (PI*x*x/2.0).cos();
	    break;
	}


	// Asymptotic power series auxiliary functions for large argument
	let x2 = x * x;
	let t = PI * x2;
	let u = 1.0 / (t * t);
	let t = 1.0 / t;
	let f = 1.0 - u * polevl(u,&f_n,9) / p1evl(u,&f_d,10); // f_n : 9 ??
	let g = t * polevl(u,&g_n,10) / p1evl(u,&g_d,11);

	let t = FRAC_PI_2 * x2;
	let c = t.cos();
	let s = t.sin();
	let t = PI * x;
	cc = 0.5 + (f * s - g * c) / t;
	ss = 0.5 - (f * c + g * s) / t;

	break;
    }

    if xxa < 0.0 {
	cc = -cc;
	ss = -ss;
    }

    (ss,cc)
}

#[cfg(test)]
mod test_data;

#[test]
fn test_fresnel() {
    let tol = R::EPSILON;
    for &[x,c1,s1] in test_data::TEST_DATA {
	let (c2,s2) = fresnl(x);
	let e = (c1 - c2).abs().max(s1 - s2).abs();
	if e > tol {
	    panic!("{x} : {c1},{s1} vs {c2},{s2} : {e:e} > {:e}",tol);
	}
    }
}
