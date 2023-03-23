#[derive(Clone, Copy)]
pub struct Scheme<O, M>
where
    O: Fn(f64, f64, f64, f64) -> (f64, f64),
    M: Fn(&Vec<f64>, &Vec<f64>, f64, f64, usize, f64, f64) -> (f64, f64),
{
    x0: f64,
    y0: f64,
    t: (f64, f64),
    h0: f64,
    one_step: O,
    multi_step: M,
}

impl<O, M> Scheme<O, M>
where
    O: Fn(f64, f64, f64, f64) -> (f64, f64),
    M: Fn(&Vec<f64>, &Vec<f64>, f64, f64, usize, f64, f64) -> (f64, f64),
{
    pub fn new(
        x0: f64,
        y0: f64,
        t: (f64, f64),
        h0: f64,
        one_step: O,
        multi_step: M,
    ) -> Scheme<O, M> {
        Scheme {
            x0,
            y0,
            t,
            h0,
            one_step,
            multi_step,
        }
    }

    fn optimal_h(&self, eps: f64) -> f64 {
        let mut h = self.h0;
        let mut t = self.t.0;
        let mut h_opt = h / 2.;

        let (mut x_next, mut y_next) = (self.one_step)(self.x0, self.y0, self.t.0, h);
        t += h_opt;

        let (x_half, y_half) = (self.one_step)(self.x0, self.y0, self.t.0, h_opt);
        t += h_opt;

        let (mut x_next_half, mut y_next_half) = (self.one_step)(x_half, y_half, t, h_opt);

        while (x_next_half - x_next).abs() + (y_next_half - y_next).abs() > eps {
            t = self.t.0;
            h = h_opt;
            h_opt /= 2.;

            t += h_opt;
            (x_next, y_next) = (self.one_step)(self.x0, self.y0, t, h);

            t += h_opt;
            let (x_half, y_half) = (self.one_step)(self.x0, self.y0, self.t.0, h_opt);

            t += h_opt;
            (x_next_half, y_next_half) = (self.one_step)(x_half, y_half, t, h_opt);
        }

        h_opt
    }

    pub fn iterate(&self, tarr: &mut Vec<f64>, eps: f64, count: usize) -> (Vec<f64>, Vec<f64>) {
        let (mut x, mut y) = (vec![self.x0], vec![self.y0]);
        let mut t = self.t.0;
        let h = self.optimal_h(eps);

        for _ in 1..(count + 1) {
            t += h;
            let (x_res, y_res) = (self.one_step)(*x.last().unwrap(), *y.last().unwrap(), t, h);
            tarr.push(t);
            x.push(x_res);
            y.push(y_res);
        }

        let mut i = count;
        while t <= self.t.1 {
            let mut x_iter = *x.last().unwrap();
            let mut y_iter = *y.last().unwrap();

            let (mut x_next, mut y_next) = (self.multi_step)(&x, &y, x_iter, y_iter, i, t, h);

            while (x_next - x_iter).abs() + (y_next - y_iter).abs() > eps {
                (x_iter, y_iter) = (x_next, y_next);

                (x_next, y_next) = (self.multi_step)(&x, &y, x_iter, y_iter, i, t, h);
            }

            x.push(x_next);
            y.push(y_next);
            t += h;
            tarr.push(t);
            i += 1;
        }

        (x, y)
    }
}
