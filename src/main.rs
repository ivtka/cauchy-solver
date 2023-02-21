mod derive;

use derive::find_optimal_h;

fn main() {
    let eps = 1e-2;
    let h0 = 1e-3;
    let t0 = 0.0;

    let x0 = 0.0;
    let y0 = 1.0;

    let mut _t = vec![t0];

    let _real_x = |t: f64| t.sin() + 1.0;

    let _real_y = |t: f64| t.powi(2);

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

    let _optimal_h = find_optimal_h(x0, y0, t0, h0, eps, yrk21);
}
