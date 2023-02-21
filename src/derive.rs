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
