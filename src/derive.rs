/// Find Optimal `H` using RK4 algorithm.
///
/// # Arguments
///
/// * `f` - One step difference scheme
///
/// # Examples
///
/// ```
/// let optimal_h = find_optimal_h(
///     x0,
///     y0,
///     Some(t0),
///     Some(h0),
///     Some(eps),
///     |h0: f64, x0: f64, y0: f64, t0: f64| {
///         let dx =
///             |t: f64, x: f64, y: f64| 2.0 * x - y * t.powi(2) - 2.0 * (t.sin() + 1.0) * t.cos();
///
///         let dy =
///             |t: f64, x: f64, y: f64| x + 2.0 * y - t.sin() - 2.0 * t.powi(2) + 2.0 * t - 1.0;
///
///         let k1 = dx(t0, x0, y0);
///         let k2 = dx(t0 + h0, x0 + h0 * k1, y0);
///
///         let l1 = dy(t0, x0, y0);
///         let l2 = dy(t0 + h0, x0, y0 + h0 * l1);
///
///         let x_next = x0 + h0 / 2.0 * (k1 + k2);
///         let y_next = y0 + h0 / 2.0 * (l1 + l2);
///
///         (x_next, y_next)
///     },
/// );
/// ```
pub fn find_optimal_h<F>(x0: f64, y0: f64, t0: f64, h0: f64, epsilon: f64, f: F) -> f64
where
    F: Fn(f64, f64, f64, f64) -> (f64, f64),
{
    let mut optimal_h = h0 / 2.0;
    let mut t = t0;

    let (mut x_next, mut y_next) = f(h0, x0, y0, t);

    t += optimal_h;

    let (x_half, y_half) = f(optimal_h, x0, y0, t);

    t += optimal_h;

    let (mut x_half_next, mut y_half_next) = f(optimal_h, x_half, y_half, t);

    while (x_half_next - x_next).abs() + (y_half_next - y_next).abs() > epsilon {
        t = t0;
        let h = optimal_h;
        optimal_h /= 2.0;

        (x_next, y_next) = f(h, x0, y0, t);

        t += optimal_h;

        let (x_half, y_half) = f(optimal_h, x0, y0, t);

        (x_half_next, y_half_next) = f(optimal_h, x_half, y_half, t);
    }

    optimal_h
}

pub struct Scheme(pub Vec<f64>, pub Vec<f64>);

pub fn newton<F>(scheme: &mut Scheme, h: f64, eps: f64, f: F)
where
    F: Fn(&Vec<f64>, usize, f64, f64, f64, f64) -> f64,
{
    let mut t = h * 3.0;
    let (mut x_s, mut y_s) = (scheme.0[2], scheme.1[2]);

    let mut matrix = [[0.0; 2]; 2];
    let mut matrix_inv = [[0.0; 2]; 2];

    let mut i = 4;

    let (mut x_next, mut y_next);

    while (t - 1.0).abs() > eps {
        // println!("t < 1.0 == {}, {t}", t < 1.0);
        loop {
            let a = 1.0 - 4.0 / 3.0 * h;
            let b = -(2.0 / 3.0) * h;
            let c = -(2.0 / 3.0) * h;
            let d = 1.0 + 2.0 / 3.0 * h;

            matrix[0][0] = a;
            matrix[0][1] = b;
            matrix[1][0] = c;
            matrix[1][1] = d;

            let det = find_determinator(&matrix);

            matrix_inv[0][0] = matrix[1][1] * (1.0 / det);
            matrix_inv[1][1] = matrix[0][0] * (1.0 / det);
            matrix_inv[0][1] = -1.0 * matrix[0][1] * (1.0 / det);
            matrix_inv[1][0] = -1.0 * matrix[1][0] * (1.0 / det);

            let f_x = x_s - f(&scheme.0, i, h, t, x_s, y_s);
            let f_y = x_s - f(&scheme.1, i, h, t, x_s, y_s);

            x_next = x_s - (matrix_inv[0][0]) * f_x + matrix_inv[0][1] * f_y;
            y_next = x_s - (matrix_inv[1][0]) * f_x + matrix_inv[0][1] * f_y;

            if (x_next - x_s).abs() + (y_next - y_s).abs() < eps {
                break;
            } else {
                (x_s, y_s) = (x_next, y_next);
            }
        }

        (scheme.0[i + 1], scheme.1[i + 1]) = (x_next, y_next);

        (x_s, y_s) = (x_next, y_next);

        t += h;
        i += 1;
    }
}

fn find_determinator(matrix: &[[f64; 2]; 2]) -> f64 {
    matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1]
}
