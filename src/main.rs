mod derive;

use derive::*;

fn main() {
    let eps = 1e-2;
    let h0 = 1e-3;
    let mut t0 = 0.0;

    let x0 = 0.0;
    let y0 = 1.0;

    let dx = |t: f64, x: f64, y: f64| 2.0 * x - y * t.powi(2) - 2.0 * (t.sin() + 1.0) * t.cos();

    let dy = |t: f64, x: f64, y: f64| x + 2.0 * y - t.sin() - 2.0 * t.powi(2) + 2.0 * t - 1.0;

    let yrk21 = |h0: f64, x0: f64, y0: f64, t0: f64| {
        let k1 = dx(t0, x0, y0);
        let k2 = dx(t0 + h0, x0 + h0 * k1, y0);

        let l1 = dy(t0, x0, y0);
        let l2 = dy(t0 + h0, x0, y0 + h0 * l1);

        let x_next = x0 + h0 / 2.0 * (k1 + k2);
        let y_next = y0 + h0 / 2.0 * (l1 + l2);

        (x_next, y_next)
    };

    let ng4 = |values: &Vec<f64>, i: usize, h: f64, t: f64, x_s: f64, y_s: f64| {
        48.0 / 25.0 * values[i - 4] - 36.0 / 25.0 * values[i - 3] + 16.0 / 25.0 * values[i - 2]
            - 3.0 / 25.0 * values[i - 1]
            + 12.0 / 25.0 * h * dx(t, x_s, y_s)
    };
    let optimal_h = find_optimal_h(x0, y0, t0, h0, eps, yrk21);

    let num = (1.0 / optimal_h + 1.0) as usize;

    let mut scheme = Scheme(vec![0.0; num], vec![0.0; num]);

    fill_scheme(&mut scheme, x0, y0, optimal_h, yrk21);
    newton(&mut scheme, optimal_h, eps, ng4);

    let mut i = 0;

    while (t0 - 1.0).abs() > eps {
        let (x_t, y_t) = (t0.powi(2).exp() / 2.0, (-t0.powi(2).exp() / 2.0));

        println!("({}, {}) - ({x_t}, {y_t})", scheme.0[i], scheme.1[i]);

        t0 += optimal_h;
        i += 1;
    }
}

fn fill_scheme<F>(scheme: &mut Scheme, x0: f64, y0: f64, h: f64, f: F)
where
    F: Fn(f64, f64, f64, f64) -> (f64, f64),
{
    let mut t = 0.0;

    let (mut x0, mut y0) = (x0, y0);

    (scheme.0[0], scheme.1[0]) = (x0, y0);

    for i in 1..3 {
        let (x_next, y_next) = f(h, x0, y0, t);

        t += h;
        (scheme.0[i], scheme.1[i]) = (x_next, y_next);
        (x0, y0) = (x_next, y_next);
    }
}
