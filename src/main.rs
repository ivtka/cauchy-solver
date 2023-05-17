mod derive;

use derive::Scheme;
use std::fs::File;
use std::io::{BufWriter, Write};

fn main() {
    let real_x = |t: f64| t.sin() + 1.;
    let real_y = |t: f64| t.powi(2);

    let yrk21 = |x0: f64, y0: f64, t0: f64, h0: f64| {
        let k1 = dx(t0, x0, y0);
        let k2 = dx(t0 + h0, x0 + h0 * k1, y0 + h0 * k1);

        let l1 = dy(t0, x0, y0);
        let l2 = dy(t0 + h0, x0 + h0 * l1, y0 + h0 * l1);

        (x0 + h0 * (k1 + k2) / 2., y0 + h0 * (l1 + l2) / 2.)
    };

    let ng4 = |x: &Vec<f64>, y: &Vec<f64>, x_iter: f64, y_iter: f64, i: usize, t: f64, h: f64| {
        let x = 48. / 25. * x[i] - 36. / 25. * x[i - 1] + 16. / 25. * x[i - 2]
            - 3. / 25. * x[i - 3]
            + 12. / 25. * h * dx(t, x_iter, y_iter);
        let y = 48. / 25. * y[i] - 36. / 25. * y[i - 1] + 16. / 25. * y[i - 2]
            - 3. / 25. * y[i - 3]
            + 12. / 25. * h * dy(t, x_iter, y_iter);

        (x, y)
    };

    let (x0, y0) = (1., 0.);
    let t = (0., 1.);
    let h0 = 1e-2;

    let scheme = Scheme::new(x0, y0, t, h0, yrk21, ng4);

    let eps = 1e-2;
    let mut tarr = vec![0.];

    let (x, y) = scheme.iterate(&mut tarr, eps, 4);

    let mut j = 0;
    let file = File::create("iterate.txt").expect("unable to create file iterate.txt");
    let mut file = BufWriter::new(file);
    for i in 0..x.len() {
        if (tarr[i] - j as f64 * eps).abs() < eps {
            writeln!(
                file,
                "({:.5}, {:.5}) - ({:.5}, {:.5})",
                x[i],
                y[i],
                real_x(tarr[i]),
                real_y(tarr[i])
            )
            .expect("unable to write");
            j += 1;
        }
    }

    let mut t0 = t.0;
    let steps = ((t.1 - t.0) / h0) as usize;

    let eps = 1e-2;

    let file = File::create("newton.txt").expect("unable to create file newton.txt");
    let mut file = BufWriter::new(file);
    for i in 0..steps {
        t0 += h0;
        let (x0, y0) = yrk21(x0, y0, t0, h0);
        let (x, y) = newton(x0, y0, t0, eps, steps);
        writeln!(
            file,
            "({:.5}, {:.5}) - ({:.5}, {:.5})",
            x,
            y,
            real_x(t0),
            real_y(t0)
        )
        .expect("unable to write");
    }
}

fn dx(t: f64, x: f64, y: f64) -> f64 {
    2. * x - y + t.powi(2) - 2. * (t.sin() + 1.) + t.cos()
}

fn dy(t: f64, x: f64, y: f64) -> f64 {
    x + 2. * y - t.sin() - 2. * t.powi(2) + 2. * t - 1.
}

fn newton(x0: f64, y0: f64, t0: f64, eps: f64, max_iter: usize) -> (f64, f64) {
    let mut x = x0;
    let mut y = y0;
    let mut i = 0;

    while i < max_iter {
        let x_res = dx(t0, x, y);
        let y_res = dy(t0, x, y);

        let j11 = 2.0;
        let j12 = -1.0;
        let j21 = 1.0;
        let j22 = 2.0;

        let det = j11 * j22 - j12 * j21;

        let x_next = (-j22 * x_res + j12 * y_res) / det;
        let y_next = (j21 * x_res - j11 * y_res) / det;

        x += x_next;
        y += y_next;

        i += 1;

        if x_res.abs() + y_res.abs() < eps {
            break;
        }
    }

    (x, y)
}
