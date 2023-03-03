#[doc = r"Find Optimal `H` using RK4 algorithm.

# Arguments

* `f` - One step difference scheme
"]
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

pub enum DeriveFn {
    Dx,
    Dy,
}

impl Scheme {
    #[doc = r"Iterate method

# Arguments

* `f` - Multi step difference scheme
"]
    pub fn iterate<F>(&mut self, h: f64, eps: f64, f: F, level: usize)
    where
        F: Fn(&Vec<f64>, usize, f64, f64, f64, f64, DeriveFn) -> f64,
    {
        let mut t = h * 3.0;

        let (mut x_s, mut y_s) = (self.0[1], self.1[1]);

        let mut i = level;

        let (mut x_next, mut y_next);

        while (t - 1.0).abs() > eps {
            loop {
                x_next = f(&self.0, i, h, t, x_s, y_s, DeriveFn::Dx);
                y_next = f(&self.1, i, h, t, x_s, y_s, DeriveFn::Dy);

                if (x_next - x_s).abs() + (y_next - y_s).abs() > eps {
                    break;
                } else {
                    (x_s, y_s) = (x_next, y_next);
                }
            }

            (self.0[i + 1], self.1[i + 1]) = (x_next, y_next);

            (x_s, y_s) = (x_next, y_next);

            t += h;
            i += 1;
        }
    }

    pub fn new<F>(x0: f64, y0: f64, h: f64, f: F, count: usize, num: usize) -> Scheme
    where
        F: Fn(f64, f64, f64, f64) -> (f64, f64),
    {
        let mut t = 0.0;

        let mut scheme = Scheme(vec![0.0; num], vec![0.0; num]);

        let (mut x0, mut y0) = (x0, y0);

        (scheme.0[0], scheme.1[0]) = (x0, y0);

        for i in 1..count {
            let (x_next, y_next) = f(h, x0, y0, t);

            t += h;

            (scheme.0[i], scheme.1[i]) = (x_next, y_next);
            (x0, y0) = (x_next, y_next);
        }

        scheme
    }

    #[doc = r"Newton method

# Arguments

* `f` - Multi step difference scheme
"]
    pub fn newton<F>(&mut self, h: f64, eps: f64, f: F, level: usize)
    where
        F: Fn(&Vec<f64>, usize, f64, f64, f64, f64, DeriveFn) -> f64,
    {
        let mut t = h * 3.0;
        let (mut x_s, mut y_s) = (self.0[2], self.1[2]);

        let mut matrix = [[0.0; 2]; 2];
        let mut matrix_inv = [[0.0; 2]; 2];

        let mut i = level;

        let (mut x_next, mut y_next);

        let determinator =
            |matrix: &[[f64; 2]; 2]| matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];

        while (t - 1.0).abs() > eps {
            loop {
                let a = 1.0 - 4.0 / 3.0 * h;
                let b = -(2.0 / 3.0) * h;
                let c = -(2.0 / 3.0) * h;
                let d = 1.0 + 2.0 / 3.0 * h;

                matrix[0][0] = a;
                matrix[0][1] = b;
                matrix[1][0] = c;
                matrix[1][1] = d;

                let det = determinator(&matrix);

                matrix_inv[0][0] = matrix[1][1] * (1.0 / det);
                matrix_inv[1][1] = matrix[0][0] * (1.0 / det);
                matrix_inv[0][1] = -1.0 * matrix[0][1] * (1.0 / det);
                matrix_inv[1][0] = -1.0 * matrix[1][0] * (1.0 / det);

                let f_x = x_s - f(&self.0, i, h, t, x_s, y_s, DeriveFn::Dx);
                let f_y = y_s - f(&self.1, i, h, t, x_s, y_s, DeriveFn::Dy);

                x_next = x_s - (matrix_inv[0][0]) * f_x + matrix_inv[0][1] * f_y;
                y_next = y_s - (matrix_inv[1][0]) * f_x + matrix_inv[0][1] * f_y;

                if (x_next - x_s).abs() + (y_next - y_s).abs() > eps {
                    break;
                } else {
                    (x_s, y_s) = (x_next, y_next);
                }
            }

            (self.0[i + 1], self.1[i + 1]) = (x_next, y_next);

            (x_s, y_s) = (x_next, y_next);

            t += h;
            i += 1;
        }
    }
}