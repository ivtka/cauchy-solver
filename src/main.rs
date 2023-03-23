mod derive;

use derive::Scheme;
use std::fs::File;
use std::io::{BufWriter, Write};

fn main() {
    let dx = |t: f64, x: f64, y: f64| 2. * x - y + t.powi(2) - 2. * (t.sin() + 1.) + t.cos();
    let dy = |t: f64, x: f64, y: f64| x + 2. * y - t.sin() - 2. * t.powi(2) + 2. * t - 1.;

    let real_x = |t: f64| t.sin() + 1.;
    let real_y = |t: f64| t.powi(2);

    let yrk21 = |x0: f64, y0: f64, t0: f64, h0: f64| {
        let k1 = dx(t0, x0, y0);
        let k2 = dx(t0 + h0, x0 + h0 * k1, y0 + h0 * k1);

        let l1 = dy(t0, x0, y0);
        let l2 = dy(t0 + h0, x0 + h0 * k1, y0 + h0 * l1);

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

    let scheme = Scheme::new(1., 0., (0., 1.), 1e-4, yrk21, ng4);

    let eps = 1e-3;
    let mut tarr = vec![0.];

    let (x, y) = scheme.iterate(&mut tarr, eps, 3);

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
}
