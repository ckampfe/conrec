/*
Derivation from the fortran version of CONREC by Paul Bourke
d               ! matrix of data to contour
ilb,iub,jlb,jub ! index bounds of data matrix
x               ! data matrix column coordinates
y               ! data matrix row coordinates
nc              ! number of contour levels
z               ! contour levels in increasing order
*/
// #[allow(clippy::too_many_arguments)]
pub fn contour<F>(
    d: &[&[f64]],
    ilb: isize,
    iub: isize,
    jlb: isize,
    jub: isize,
    x: &[f64],
    y: &[f64],
    nc: isize,
    z: &[f64],
    mut conrec_line: F,
) where
    F: FnMut(f64, f64, f64, f64, f64),
{
    let mut m1: isize = 0;
    let mut m2: isize = 0;
    let mut m3: isize = 0;
    let mut case_value: isize = 0;
    let mut dmin: f64 = 0.0;
    let mut dmax: f64 = 0.0;
    let mut x1: f64 = 0.0;
    let mut x2: f64 = 0.0;
    let mut y1: f64 = 0.0;
    let mut y2: f64 = 0.0;
    let mut i: isize = 0;
    let mut j: isize = 0;
    let mut k: isize = 0;
    let mut m: isize = 0;
    let mut h: [f64; 5] = [0.0; 5];
    let mut sh: [isize; 5] = [0; 5];
    let mut xh: [f64; 5] = [0.0; 5];
    let mut yh: [f64; 5] = [0.0; 5];
    let im: [isize; 4] = [0, 1, 1, 0];
    let jm: [isize; 4] = [0, 0, 1, 1];
    let castab: [[[isize; 3]; 3]; 3] = [
        [[0, 0, 8], [0, 2, 5], [7, 6, 9]],
        [[0, 3, 4], [1, 3, 1], [4, 3, 0]],
        [[9, 6, 7], [5, 2, 0], [8, 0, 0]],
    ];
    let mut temp1: f64 = 0.;
    let mut temp2: f64 = 0.;
    j = jub - 1;
    while j >= jlb {
        i = ilb;
        while i < iub {
            temp1 = d[i as usize][j as usize].min(d[i as usize][j as usize + 1]);
            temp2 = d[i as usize + 1][j as usize].min(d[i as usize + 1][j as usize + 1]);
            dmin = temp1.min(temp2);
            temp1 = d[i as usize][j as usize].max(d[i as usize][j as usize + 1]);
            temp2 = d[i as usize + 1][j as usize].max(d[i as usize + 1][j as usize + 1]);
            dmax = temp1.max(temp2);

            if !(dmax < z[0] || dmin > z[nc as usize - 1]) {
                k = 0;
                while k < nc {
                    if !(z[k as usize] < dmin || z[k as usize] > dmax) {
                        m = 4;
                        while m >= 0 {
                            if m > 0 {
                                h[m as usize] = d[i as usize + im[m as usize - 1] as usize]
                                    [j as usize + jm[m as usize - 1] as usize]
                                    - z[k as usize];
                                xh[m as usize] = x[i as usize + im[m as usize - 1] as usize];
                                yh[m as usize] = y[j as usize + jm[m as usize - 1] as usize];
                            } else {
                                h[0] = 0.25f64 * (h[1] + h[2] + h[3] + h[4]);
                                xh[0] = 0.50 * (x[i as usize] + x[i as usize + 1]);
                                yh[0] = 0.50 * (y[j as usize] + y[j as usize + 1]);
                            }
                            if h[m as usize] > 0.0f64 {
                                sh[m as usize] = 1
                            } else if h[m as usize] < 0.0f64 {
                                sh[m as usize] = -1
                            } else {
                                sh[m as usize] = 0
                            }
                            m -= 1
                        }
                        /* k - contour */
                        /*
                           Note: at this stage the relative heights of the corners and the
                           centre are in the h array, and the corresponding coordinates are
                           in the xh and yh arrays. The centre of the box is indexed by 0
                           and the 4 corners by 1 to 4 as shown below.
                           Each triangle is then indexed by the parameter m, and the 3
                           vertices of each triangle are indexed by parameters m1,m2,and m3.
                           It is assumed that the centre of the box is always vertex 2
                           though this isimportant only when all 3 vertices lie exactly on
                           the same contour level, in which case only the side of the box
                           is drawn.
                              vertex 4 +-------------------+ vertex 3
                                       | \               / |
                                       |   \    m-3    /   |
                                       |     \       /     |
                                       |       \   /       |
                                       |  m=2    X   m=2   |       the centre is vertex 0
                                       |       /   \       |
                                       |     /       \     |
                                       |   /    m=1    \   |
                                       | /               \ |
                              vertex 1 +-------------------+ vertex 2
                        */
                        /* Scan each triangle in the box */
                        m = 1;
                        while m <= 4 {
                            m1 = m;
                            m2 = 0;
                            if m != 4 {
                                m3 = m + 1
                            } else {
                                m3 = 1
                            }
                            case_value = castab[(sh[m1 as usize] + 1) as usize]
                                [(sh[m2 as usize] + 1) as usize]
                                [(sh[m3 as usize] + 1) as usize];
                            if case_value != 0 {
                                match case_value {
                                    1 => {
                                        /* Line between vertices 1 and 2 */
                                        x1 = xh[m1 as usize];
                                        y1 = yh[m1 as usize];
                                        x2 = xh[m2 as usize];
                                        y2 = yh[m2 as usize]
                                    }
                                    2 => {
                                        /* Line between vertices 2 and 3 */
                                        x1 = xh[m2 as usize];
                                        y1 = yh[m2 as usize];
                                        x2 = xh[m3 as usize];
                                        y2 = yh[m3 as usize]
                                    }
                                    3 => {
                                        /* Line between vertices 3 and 1 */
                                        x1 = xh[m3 as usize];
                                        y1 = yh[m3 as usize];
                                        x2 = xh[m1 as usize];
                                        y2 = yh[m1 as usize]
                                    }
                                    4 => {
                                        /* Line between vertex 1 and side 2-3 */
                                        x1 = xh[m1 as usize];
                                        y1 = yh[m1 as usize];
                                        x2 = (h[m3 as usize] * xh[m2 as usize]
                                            - h[m2 as usize] * xh[m3 as usize])
                                            / (h[m3 as usize] - h[m2 as usize]);
                                        y2 = (h[m3 as usize] * yh[m2 as usize]
                                            - h[m2 as usize] * yh[m3 as usize])
                                            / (h[m3 as usize] - h[m2 as usize])
                                    }
                                    5 => {
                                        /* Line between vertex 2 and side 3-1 */
                                        x1 = xh[m2 as usize];
                                        y1 = yh[m2 as usize];
                                        x2 = (h[m1 as usize] * xh[m3 as usize]
                                            - h[m3 as usize] * xh[m1 as usize])
                                            / (h[m1 as usize] - h[m3 as usize]);
                                        y2 = (h[m1 as usize] * yh[m3 as usize]
                                            - h[m3 as usize] * yh[m1 as usize])
                                            / (h[m1 as usize] - h[m3 as usize])
                                    }
                                    6 => {
                                        /* Line between vertex 3 and side 1-2 */
                                        x1 = xh[m3 as usize];
                                        y1 = yh[m3 as usize];
                                        x2 = (h[m2 as usize] * xh[m1 as usize]
                                            - h[m1 as usize] * xh[m2 as usize])
                                            / (h[m2 as usize] - h[m1 as usize]);
                                        y2 = (h[m2 as usize] * yh[m1 as usize]
                                            - h[m1 as usize] * yh[m2 as usize])
                                            / (h[m2 as usize] - h[m1 as usize])
                                    }
                                    7 => {
                                        /* Line between sides 1-2 and 2-3 */
                                        x1 = (h[m2 as usize] * xh[m1 as usize]
                                            - h[m1 as usize] * xh[m2 as usize])
                                            / (h[m2 as usize] - h[m1 as usize]);
                                        y1 = (h[m2 as usize] * yh[m1 as usize]
                                            - h[m1 as usize] * yh[m2 as usize])
                                            / (h[m2 as usize] - h[m1 as usize]);
                                        x2 = (h[m3 as usize] * xh[m2 as usize]
                                            - h[m2 as usize] * xh[m3 as usize])
                                            / (h[m3 as usize] - h[m2 as usize]);
                                        y2 = (h[m3 as usize] * yh[m2 as usize]
                                            - h[m2 as usize] * yh[m3 as usize])
                                            / (h[m3 as usize] - h[m2 as usize])
                                    }
                                    8 => {
                                        /* Line between sides 2-3 and 3-1 */
                                        x1 = (h[m3 as usize] * xh[m2 as usize]
                                            - h[m2 as usize] * xh[m3 as usize])
                                            / (h[m3 as usize] - h[m2 as usize]);
                                        y1 = (h[m3 as usize] * yh[m2 as usize]
                                            - h[m2 as usize] * yh[m3 as usize])
                                            / (h[m3 as usize] - h[m2 as usize]);
                                        x2 = (h[m1 as usize] * xh[m3 as usize]
                                            - h[m3 as usize] * xh[m1 as usize])
                                            / (h[m1 as usize] - h[m3 as usize]);
                                        y2 = (h[m1 as usize] * yh[m3 as usize]
                                            - h[m3 as usize] * yh[m1 as usize])
                                            / (h[m1 as usize] - h[m3 as usize])
                                    }
                                    9 => {
                                        /* Line between sides 3-1 and 1-2 */
                                        x1 = (h[m1 as usize] * xh[m3 as usize]
                                            - h[m3 as usize] * xh[m1 as usize])
                                            / (h[m1 as usize] - h[m3 as usize]);
                                        y1 = (h[m1 as usize] * yh[m3 as usize]
                                            - h[m3 as usize] * yh[m1 as usize])
                                            / (h[m1 as usize] - h[m3 as usize]);
                                        x2 = (h[m2 as usize] * xh[m1 as usize]
                                            - h[m1 as usize] * xh[m2 as usize])
                                            / (h[m2 as usize] - h[m1 as usize]);
                                        y2 = (h[m2 as usize] * yh[m1 as usize]
                                            - h[m1 as usize] * yh[m2 as usize])
                                            / (h[m2 as usize] - h[m1 as usize])
                                    }
                                    _ => {}
                                }
                                /* Finally draw the line */
                                conrec_line(x1, y1, x2, y2, z[k as usize]);
                            }
                            m += 1
                        }
                    }
                    k += 1
                    /* m */
                }
            }
            i += 1
        }
        j -= 1
        /* i */
    }
    /* j */
}

#[cfg(test)]
mod tests {
    // use super::*;

    // #[test]
    // fn it_works() {
    //     contour(
    //         &[
    //             &[0., 0., 0., 0., 0.],
    //             &[0., 3., 3., 3., 0.],
    //             &[0., 3., 5., 3., 0.],
    //             &[0., 3., 3., 3., 0.],
    //             &[0., 0., 0., 0., 0.],
    //         ],
    //         0,
    //         4,
    //         0,
    //         4,
    //         &[0., 1., 2., 3., 4.],
    //         &[0., 1., 2., 3., 4.],
    //         1,
    //         &[3.0],
    //         |x1, y1, x2, y2, level| println!("({x1}, {y1}) - ({x2}, {y2}), {level}"),
    //     );
    // }
}
