/*******************************************************************************

        These routines are to calculate the transpose of interpolation matrices

        Created by Se Young Chun, Jan 19, 2007, the University of Michigan

*******************************************************************************/

#include	<math.h>
#include	<stddef.h>
#include	<stdio.h>
#include	<stdlib.h>

#include	"interpol2Dtran.h"
 
extern double	InterpolatedValueTranMirr (
                        float   *Bout,      
                        long    Width,      
                        long    Height,    
                        float   Bcoeff,   
                        double  x,        
                        double  y,        
                        long    SplineDegree 
		)
{ /* begin InterpolatedValueTranMirr */

	float	*p;
	double	xWeight[10], yWeight[10];
	double	w, w2, w4, t, t0, t1;
	long	xIndex[10], yIndex[10];
	long	Width2 = 2L * Width - 2L, Height2 = 2L * Height - 2L;
	long	i, j, k;

	/* compute the interpolation indexes */
	if (SplineDegree & 1L) {
		i = (long)floor(x) - SplineDegree / 2L;
		j = (long)floor(y) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
		}
	}
	else {
		i = (long)floor(x + 0.5) - SplineDegree / 2L;
		j = (long)floor(y + 0.5) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
		}
	}

	/* compute the interpolation weights */
	switch (SplineDegree) {
		case 2L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[1] = 3.0 / 4.0 - w * w;
			xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
			xWeight[0] = 1.0 - xWeight[1] - xWeight[2];
			/* y */
			w = y - (double)yIndex[1];
			yWeight[1] = 3.0 / 4.0 - w * w;
			yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
			yWeight[0] = 1.0 - yWeight[1] - yWeight[2];
			break;
		case 3L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[3] = (1.0 / 6.0) * w * w * w;
			xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- xWeight[3];
			xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
			xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];
			/* y */
			w = y - (double)yIndex[1];
			yWeight[3] = (1.0 / 6.0) * w * w * w;
			yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- yWeight[3];
			yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
			yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];
			break;
		case 4L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= (1.0 / 24.0) * xWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			xWeight[1] = t1 + t0;
			xWeight[3] = t1 - t0;
			xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
			xWeight[2] = 1.0 - xWeight[0] - xWeight[1] 
				- xWeight[3] - xWeight[4];
			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= (1.0 / 24.0) * yWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			yWeight[1] = t1 + t0;
			yWeight[3] = t1 - t0;
			yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
			yWeight[2] = 1.0 - yWeight[0] - yWeight[1] 
				- yWeight[3] - yWeight[4];
			break;
		case 5L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- xWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			xWeight[2] = t0 + t1;
			xWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			xWeight[1] = t0 + t1;
			xWeight[4] = t0 - t1;
			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- yWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			yWeight[2] = t0 + t1;
			yWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			yWeight[1] = t0 + t1;
			yWeight[4] = t0 - t1;
			break;
		case 6L:
			/* x */
			w = x - (double)xIndex[3];
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0] * xWeight[0];
			xWeight[0] *= xWeight[0] / 720.0;
			xWeight[1] = (361.0 / 192.0 - w * (59.0 / 8.0 + w
				* (-185.0 / 16.0 + w * (25.0 / 3.0 
				+ w * (-5.0 / 2.0 + w)
				* (1.0 / 2.0 + w))))) / 120.0;
			xWeight[2] = (10543.0 / 960.0 + w * (-289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (43.0 / 6.0 
				+ w * (-17.0 / 4.0 + w
				* (-1.0 + w)))))) / 48.0;
			w2 = w * w;
			xWeight[3] = (5887.0 / 320.0 - w2 * (231.0 / 16.0 - w2
				* (21.0 / 4.0 - w2))) / 36.0;
			xWeight[4] = (10543.0 / 960.0 + w * (289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (-43.0 / 6.0 
				+ w * (-17.0 / 4.0 + w
				* (1.0 + w)))))) / 48.0;
			xWeight[6] = 1.0 / 2.0 + w;
			xWeight[6] *= xWeight[6] * xWeight[6];
			xWeight[6] *= xWeight[6] / 720.0;
			xWeight[5] = 1.0 - xWeight[0] - xWeight[1] 
				- xWeight[2] - xWeight[3]
				- xWeight[4] - xWeight[6];
			/* y */
			w = y - (double)yIndex[3];
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0] * yWeight[0];
			yWeight[0] *= yWeight[0] / 720.0;
			yWeight[1] = (361.0 / 192.0 - w * (59.0 / 8.0 + w
				* (-185.0 / 16.0 + w * (25.0 / 3.0 
				+ w * (-5.0 / 2.0 + w)
				* (1.0 / 2.0 + w))))) / 120.0;
			yWeight[2] = (10543.0 / 960.0 + w * (-289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (43.0 / 6.0 
				+ w * (-17.0 / 4.0 + w
				* (-1.0 + w)))))) / 48.0;
			w2 = w * w;
			yWeight[3] = (5887.0 / 320.0 - w2 * (231.0 / 16.0 - w2
				* (21.0 / 4.0 - w2))) / 36.0;
			yWeight[4] = (10543.0 / 960.0 + w * (289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (-43.0 / 6.0 
				+ w * (-17.0 / 4.0 + w * (1.0 + w)))))) / 48.0;
			yWeight[6] = 1.0 / 2.0 + w;
			yWeight[6] *= yWeight[6] * yWeight[6];
			yWeight[6] *= yWeight[6] / 720.0;
			yWeight[5] = 1.0 - yWeight[0] - yWeight[1] - yWeight[2] 
				- yWeight[3] - yWeight[4] - yWeight[6];
			break;
		case 7L:
			/* x */
			w = x - (double)xIndex[3];
			xWeight[0] = 1.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0] * xWeight[0];
			xWeight[0] *= (1.0 - w) / 5040.0;
			w2 = w * w;
			xWeight[1] = (120.0 / 7.0 + w * (-56.0 + w * (72.0 + w
				* (-40.0 + w2*(12.0 + w*(-6.0 + w)))))) / 720.0;
			xWeight[2] = (397.0/7.0 - w*(245.0/3.0 + w * (-15.0 + w
				* (-95.0 / 3.0 + w * (15.0 + w * (5.0 + w
				* (-5.0 + w))))))) / 240.0;
			xWeight[3] = (2416.0 / 35.0 + w2*(-48.0 + w2*(16.0 + w2
				* (-4.0 + w)))) / 144.0;
			xWeight[4] = (1191.0 / 35.0 - w * (-49.0 + w * (-9.0 + w
				* (19.0 + w*(-3.0 + w)*(-3.0 + w2))))) / 144.0;
			xWeight[5] = (40.0 / 7.0 + w*(56.0 / 3.0 + w * (24.0 + w
				*(40.0/3.0 + w2*(-4.0 + w*(-2.0 + w))))))/240.0;
			xWeight[7] = w2;
			xWeight[7] *= xWeight[7] * xWeight[7];
			xWeight[7] *= w / 5040.0;
			xWeight[6] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] 
				- xWeight[3]
				- xWeight[4] - xWeight[5] - xWeight[7];
			/* y */
			w = y - (double)yIndex[3];
			yWeight[0] = 1.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0] * yWeight[0];
			yWeight[0] *= (1.0 - w) / 5040.0;
			w2 = w * w;
			yWeight[1] = (120.0 / 7.0 + w * (-56.0 + w * (72.0 + w
				*(-40.0 + w2*(12.0 + w*(-6.0 + w)))))) / 720.0;
			yWeight[2] = (397.0 / 7.0 - w*(245.0/3.0 + w*(-15.0 + w
				* (-95.0 / 3.0 + w * (15.0 + w * (5.0 + w
				* (-5.0 + w))))))) / 240.0;
			yWeight[3] = (2416.0/35.0 + w2*(-48.0 + w2 * (16.0 + w2
				* (-4.0 + w)))) / 144.0;
			yWeight[4] = (1191.0 / 35.0 - w * (-49.0 + w * (-9.0 + w
				* (19.0 + w*(-3.0 + w)*(-3.0 + w2))))) / 144.0;
			yWeight[5] = (40.0/7.0 + w * (56.0 / 3.0 + w * (24.0 + w
				* (40.0 / 3.0 + w2 * (-4.0 + w * (-2.0 
				+ w)))))) / 240.0;
			yWeight[7] = w2;
			yWeight[7] *= yWeight[7] * yWeight[7];
			yWeight[7] *= w / 5040.0;
			yWeight[6] = 1.0 - yWeight[0] - yWeight[1] - yWeight[2] 
				- yWeight[3]
				- yWeight[4] - yWeight[5] - yWeight[7];
			break;
		case 8L:
			/* x */
			w = x - (double)xIndex[4];
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0] / 40320.0;
			w2 = w * w;
			xWeight[1] = (39.0/16.0 - w*(6.0 + w*(-9.0 / 2.0 + w2)))
				* (21.0 / 16.0 + w * (-15.0/4.0 + w*(9.0/2.0 + w
				* (-3.0 + w)))) / 5040.0;
			xWeight[2] = (82903.0 / 1792.0 + w * (-4177.0 / 32.0 + w
				* (2275.0 / 16.0 + w * (-487.0 / 8.0 
				+ w * (-85.0 / 8.0 + w
				* (41.0 / 2.0 + w * (-5.0 + w * (-2.0 
				+ w)))))))) / 1440.0;
			xWeight[3] = (310661.0/1792.0 - w * (14219.0 / 64.0 + w
				* (-199.0 / 8.0 + w * (-1327.0 / 16.0 
				+ w * (245.0 / 8.0 + w
				* (53.0 / 4.0 + w * (-8.0 + w * (-1.0 
				+ w)))))))) / 720.0;
			xWeight[4] = (2337507.0/8960.0 + w2*(-2601.0 / 16.0 + w2
				* (387.0 / 8.0 + w2 * (-9.0 + w2)))) / 576.0;
			xWeight[5] = (310661.0 / 1792.0 - w*(-14219.0 / 64.0 + w
				* (-199.0 / 8.0 + w * (1327.0 / 16.0 
				+ w * (245.0 / 8.0 + w
				* (-53.0 / 4.0 + w * (-8.0 + w * (1.0 
				+ w)))))))) / 720.0;
			xWeight[7] = (39.0/16.0 - w*(-6.0 + w*(-9.0/2.0 + w2)))
				* (21.0/16.0 + w*(15.0/4.0 + w * (9.0 / 2.0 + w
				* (3.0 + w)))) / 5040.0;
			xWeight[8] = 1.0 / 2.0 + w;
			xWeight[8] *= xWeight[8];
			xWeight[8] *= xWeight[8];
			xWeight[8] *= xWeight[8] / 40320.0;
			xWeight[6] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] 
				- xWeight[3] - xWeight[4] - xWeight[5] 
				- xWeight[7] - xWeight[8];
			/* y */
			w = y - (double)yIndex[4];
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0] / 40320.0;
			w2 = w * w;
			yWeight[1] = (39.0 / 16.0 - w*(6.0 + w*(-9.0/2.0 + w2)))
				* (21.0/16.0 + w*(-15.0/4.0 + w * (9.0 / 2.0 + w
				* (-3.0 + w)))) / 5040.0;
			yWeight[2] = (82903.0 / 1792.0 + w * (-4177.0 / 32.0 + w
				* (2275.0 / 16.0 + w * (-487.0 / 8.0 
				+ w * (-85.0 / 8.0 + w
				* (41.0 / 2.0 + w * (-5.0 + w * (-2.0 
				+ w)))))))) / 1440.0;
			yWeight[3] = (310661.0/1792.0 - w * (14219.0 / 64.0 + w
				* (-199.0 / 8.0 + w * (-1327.0 / 16.0 
				+ w * (245.0 / 8.0 + w
				* (53.0 / 4.0 + w * (-8.0 + w * (-1.0 
				+ w)))))))) / 720.0;
			yWeight[4] = (2337507.0/8960.0 + w2*(-2601.0 / 16.0 + w2
				* (387.0 / 8.0 + w2 * (-9.0 + w2)))) / 576.0;
			yWeight[5] = (310661.0 / 1792.0 - w*(-14219.0 / 64.0 + w
				* (-199.0 / 8.0 + w * (1327.0 / 16.0 
				+ w * (245.0 / 8.0 + w
				* (-53.0 / 4.0 + w * (-8.0 + w * (1.0 
				+ w)))))))) / 720.0;
			yWeight[7] = (39.0/16.0 - w*(-6.0 + w*(-9.0/2.0 + w2)))
				* (21.0 / 16.0 + w * (15.0/4.0 + w*(9.0/2.0 + w
				* (3.0 + w)))) / 5040.0;
			yWeight[8] = 1.0 / 2.0 + w;
			yWeight[8] *= yWeight[8];
			yWeight[8] *= yWeight[8];
			yWeight[8] *= yWeight[8] / 40320.0;
			yWeight[6] = 1.0 - yWeight[0] - yWeight[1] - yWeight[2] 
				- yWeight[3] - yWeight[4] - yWeight[5] 
				- yWeight[7] - yWeight[8];
			break;
		case 9L:
			/* x */
			w = x - (double)xIndex[4];
			xWeight[0] = 1.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0] * (1.0 - w) / 362880.0;
			xWeight[1] = (502.0 / 9.0 + w * (-246.0 + w * (472.0 + w
				* (-504.0 + w * (308.0 + w * (-84.0 
				+ w * (-56.0 / 3.0 + w
				* (24.0 + w * (-8.0 + w))))))))) / 40320.0;
			xWeight[2] = (3652.0 / 9.0 - w * (2023.0 / 2.0 
				+ w * (-952.0 + w
				* (938.0 / 3.0 + w * (112.0 + w * (-119.0 
				+ w * (56.0 / 3.0 + w
				* (14.0 + w * (-7.0 + w))))))))) / 10080.0;
			xWeight[3] = (44117.0 / 42.0 + w * (-2427.0 / 2.0 
				+ w * (66.0 + w
				* (434.0 + w * (-129.0 + w * (-69.0 
				+ w * (34.0 + w * (6.0 + w
				* (-6.0 + w))))))))) / 4320.0;
			w2 = w * w;
			xWeight[4] = (78095.0 / 63.0 - w2 * (700.0 
				+ w2 * (-190.0 + w2
				* (100.0 / 3.0 + w2 * (-5.0 + w))))) / 2880.0;
			xWeight[5] = (44117.0/63.0 + w * (809.0 + w * (44.0 + w
				* (-868.0 / 3.0 + w * (-86.0 + w * (46.0 
				+ w * (68.0 / 3.0 + w
				* (-4.0 + w * (-4.0 + w))))))))) / 2880.0;
			xWeight[6] = (3652.0 / 21.0 - w * (-867.0 / 2.0 
				+ w * (-408.0 + w
				* (-134.0 + w * (48.0 + w * (51.0 
				+ w * (-4.0 + w) * (-1.0 + w)
				* (2.0 + w))))))) / 4320.0;
			xWeight[7] = (251.0 / 18.0 + w * (123.0 / 2.0 
				+ w * (118.0 + w
				*(126.0 + w*(77.0 + w*(21.0 + w*(-14.0 / 3.0 + w
				* (-6.0 + w * (-2.0 + w))))))))) / 10080.0;
			xWeight[9] = w2 * w2;
			xWeight[9] *= xWeight[9] * w / 362880.0;
			xWeight[8] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] 
				- xWeight[3] - xWeight[4] - xWeight[5] 
				- xWeight[6] - xWeight[7] - xWeight[9];
			/* y */
			w = y - (double)yIndex[4];
			yWeight[0] = 1.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0] * (1.0 - w) / 362880.0;
			yWeight[1] = (502.0 / 9.0 + w * (-246.0 + w * (472.0 + w
				* (-504.0 + w * (308.0 + w * (-84.0 
				+ w * (-56.0 / 3.0 + w
				* (24.0 + w * (-8.0 + w))))))))) / 40320.0;
			yWeight[2] = (3652.0 / 9.0 - w * (2023.0 / 2.0 
				+ w * (-952.0 + w
				* (938.0 / 3.0 + w * (112.0 + w * (-119.0 
				+ w * (56.0 / 3.0 + w
				* (14.0 + w * (-7.0 + w))))))))) / 10080.0;
			yWeight[3] = (44117.0 / 42.0 + w * (-2427.0 / 2.0 
				+ w * (66.0 + w
				* (434.0 + w * (-129.0 + w * (-69.0 
				+ w * (34.0 + w * (6.0 + w
				* (-6.0 + w))))))))) / 4320.0;
			w2 = w * w;
			yWeight[4] = (78095.0/63.0 - w2*(700.0 + w2*(-190.0 + w2
				* (100.0 / 3.0 + w2 * (-5.0 + w))))) / 2880.0;
			yWeight[5] = (44117.0 / 63.0 + w * (809.0 + w*(44.0 + w
				* (-868.0 / 3.0 + w * (-86.0 + w * (46.0 
				+ w * (68.0 / 3.0 + w
				* (-4.0 + w * (-4.0 + w))))))))) / 2880.0;
			yWeight[6] = (3652.0 / 21.0 - w * (-867.0 / 2.0 
				+ w * (-408.0 + w
				* (-134.0 + w * (48.0 + w * (51.0 
				+ w * (-4.0 + w) * (-1.0 + w)
				* (2.0 + w))))))) / 4320.0;
			yWeight[7] = (251.0/18.0 + w*(123.0/2.0 + w * (118.0 + w
				* (126.0 + w*(77.0 + w*(21.0 + w*(-14.0/3.0 + w
				* (-6.0 + w * (-2.0 + w))))))))) / 10080.0;
			yWeight[9] = w2 * w2;
			yWeight[9] *= yWeight[9] * w / 362880.0;
			yWeight[8] = 1.0 - yWeight[0] - yWeight[1] - yWeight[2] 
				- yWeight[3] - yWeight[4] - yWeight[5] 
				- yWeight[6] - yWeight[7] - yWeight[9];
			break;
		default:
			printf("Invalid spline degree\n");
			return(0.0);
	}

	/* apply the mirror boundary conditions */
	for (k = 0L; k <= SplineDegree; k++) {
		xIndex[k] = (Width == 1L) ? (0L) : ((xIndex[k] < 0L) ?
			(-xIndex[k] - Width2 * ((-xIndex[k]) / Width2))
			: (xIndex[k] - Width2 * (xIndex[k] / Width2)));
		if (Width <= xIndex[k]) {
			xIndex[k] = Width2 - xIndex[k];
		}
		yIndex[k] = (Height == 1L) ? (0L) : ((yIndex[k] < 0L) ?
			(-yIndex[k] - Height2 * ((-yIndex[k]) / Height2))
			: (yIndex[k] - Height2 * (yIndex[k] / Height2)));
		if (Height <= yIndex[k]) {
			yIndex[k] = Height2 - yIndex[k];
		}
	}

	/* perform interpolation */
	for (j = 0L; j <= SplineDegree; j++) {
		p = Bout + (ptrdiff_t)(yIndex[j] * Width);

		for (i = 0L; i <= SplineDegree; i++) {
			p[xIndex[i]] += yWeight[j] * xWeight[i] * Bcoeff;
		}
	}

	return 1;

} /* end InterpolatedValueTranMirr */

extern double	InterpolatedGradXTranMirr (
                        float   *Bout,     
                        long    Width,   
                        long    Height,  
                        float   Bcoeff,     
                        double  x,       
                        double  y,       
                        long    SplineDegree  
		)
{ /* begin InterpolatedGradXTranMirr */

	float	*p;
	double	gradxWeight[10], yWeight[10];
	double	w, w2, w4, t, t0, t1;
	long	xIndex[10], yIndex[10];
	long	Width2 = 2L * Width - 2L, Height2 = 2L * Height - 2L;
	long	i, j, k;

	/* compute the interpolation indexes */
	if (SplineDegree & 1L) {
		i = (long)floor(x) - SplineDegree / 2L;
		j = (long)floor(y) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
		}
	}
	else {
		i = (long)floor(x + 0.5) - SplineDegree / 2L;
		j = (long)floor(y + 0.5) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
		}
	}

	/* compute the interpolation weights */
	switch (SplineDegree) {
		case 2L:
			/* gradx */
                        w = x - (double)xIndex[1];
                        gradxWeight[1] = - 2.0 * w;
                        gradxWeight[2] = (1.0 / 2.0) * (1.0 - gradxWeight[1]);
                        gradxWeight[0] = - gradxWeight[1] - gradxWeight[2];
			/* y */
			w = y - (double)yIndex[1];
			yWeight[1] = 3.0 / 4.0 - w * w;
			yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
			yWeight[0] = 1.0 - yWeight[1] - yWeight[2];
			break;
		case 3L:
 			/* gradx */
                        w = x - (double)xIndex[1];
                        gradxWeight[3] = (1.0 / 2.0) * w * w;
                        gradxWeight[0] = w - 1.0 / 2.0 - gradxWeight[3];
                        gradxWeight[2] = 1.0 + gradxWeight[0] 
				- 2.0 * gradxWeight[3];
                        gradxWeight[1] = - gradxWeight[0] - gradxWeight[2] 					- gradxWeight[3];
			/* y */
			w = y - (double)yIndex[1];
			yWeight[3] = (1.0 / 6.0) * w * w * w;
			yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- yWeight[3];
			yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
			yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];
			break;
		case 4L:
 			/* gradx */
                        w = x - (double)xIndex[2];
                        w2 = w * w;
                        gradxWeight[0] = - (1.0 / 8.0) * (w2 - w + 1.0 / 4.0);
                        t = (1.0 / 2.0) * w;
                        t0 = t * w - 11.0 / 24.0;
                        t1 = t - (2.0 / 3.0) * w2 * w;
                        gradxWeight[1] = t1 + t0;
                        gradxWeight[3] = t1 - t0;
                        gradxWeight[4] = gradxWeight[0] + t0 + (1.0 / 2.0);
                        gradxWeight[2] = - gradxWeight[0] - gradxWeight[1] 
				- gradxWeight[3] - gradxWeight[4];
			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= (1.0 / 24.0) * yWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			yWeight[1] = t1 + t0;
			yWeight[3] = t1 - t0;
			yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
			yWeight[2] = 1.0 - yWeight[0] - yWeight[1] 
				- yWeight[3] - yWeight[4];
			break;
		case 5L:
  			/* gradx */
                        w = x - (double)xIndex[2];
                        w2 = w * w;
                        gradxWeight[5] = (1.0 / 24.0) * w2 * w2;
                        w2 -= w;
                        w4 = w2 * w2;
                        w -= 1.0 / 2.0;
                        t = w2 * (w2 - 3.0);
                        gradxWeight[0] = (1.0 / 12.0) * w * (1.0 + 2.0 * w2) 
				- gradxWeight[5];
                        t0 = (1.0 / 12.0) * (2.0 * w2 - 5.0) * w;
                        t1 = (-1.0 / 12.0) * (w4 - 3.0 * w2 + 4.0) 
				- (1.0 / 6.0) * w * w * (2.0 * w2 - 3.0);
                        gradxWeight[2] = t0 + t1;
                        gradxWeight[3] = t0 - t1;
                        t0 = (1.0 / 8.0) * (3.0 - 2.0 * w2) * w;
                        t1 = (1.0 / 24.0) * (w4 - w2 - 5.0 
				+ 2.0 * w * w * (2.0 * w2 - 1));
                        gradxWeight[1] = t0 + t1;
                        gradxWeight[4] = t0 - t1;
			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- yWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			yWeight[2] = t0 + t1;
			yWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			yWeight[1] = t0 + t1;
			yWeight[4] = t0 - t1;
			break;
		default:
			printf("Invalid spline degree\n");
			return(0.0);
	}

	/* apply the mirror boundary conditions */
	for (k = 0L; k <= SplineDegree; k++) {
		xIndex[k] = (Width == 1L) ? (0L) : ((xIndex[k] < 0L) ?
			(-xIndex[k] - Width2 * ((-xIndex[k]) / Width2))
			: (xIndex[k] - Width2 * (xIndex[k] / Width2)));
		if (Width <= xIndex[k]) {
			xIndex[k] = Width2 - xIndex[k];
		}
		yIndex[k] = (Height == 1L) ? (0L) : ((yIndex[k] < 0L) ?
			(-yIndex[k] - Height2 * ((-yIndex[k]) / Height2))
			: (yIndex[k] - Height2 * (yIndex[k] / Height2)));
		if (Height <= yIndex[k]) {
			yIndex[k] = Height2 - yIndex[k];
		}
	}

	/* perform interpolation */
	for (j = 0L; j <= SplineDegree; j++) {
		p = Bout + (ptrdiff_t)(yIndex[j] * Width);

		for (i = 0L; i <= SplineDegree; i++) {
			p[xIndex[i]] += yWeight[j] * gradxWeight[i] * Bcoeff;
		}
	}

	return 1;

} /* end InterpolatedGradXTranMirr */

extern double	InterpolatedGradYTranMirr (
                        float   *Bout,     
                        long    Width,    
                        long    Height,   
                        float   Bcoeff,  
                        double  x,    
                        double  y,    
                        long    SplineDegree   
		)
{ /* begin InterpolatedGradYTranMirr */

	float	*p;
	double	xWeight[10], gradyWeight[10];
	double	w, w2, w4, t, t0, t1;
	long	xIndex[10], yIndex[10];
	long	Width2 = 2L * Width - 2L, Height2 = 2L * Height - 2L;
	long	i, j, k;

	/* compute the interpolation indexes */
	if (SplineDegree & 1L) {
		i = (long)floor(x) - SplineDegree / 2L;
		j = (long)floor(y) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
		}
	}
	else {
		i = (long)floor(x + 0.5) - SplineDegree / 2L;
		j = (long)floor(y + 0.5) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
		}
	}

	/* compute the interpolation weights */
	switch (SplineDegree) {
		case 2L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[1] = 3.0 / 4.0 - w * w;
			xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
			xWeight[0] = 1.0 - xWeight[1] - xWeight[2];
 			/* grady */
                        w = y - (double)yIndex[1];
                        gradyWeight[1] = - 2.0 * w;
                        gradyWeight[2] = (1.0 / 2.0) * (1.0 - gradyWeight[1]);
                        gradyWeight[0] = - gradyWeight[1] - gradyWeight[2];
			break;
		case 3L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[3] = (1.0 / 6.0) * w * w * w;
			xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- xWeight[3];
			xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
			xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];
			/* grady */
                        w = y - (double)yIndex[1];
                        gradyWeight[3] = (1.0 / 2.0) * w * w;
                        gradyWeight[0] = w - 1.0 / 2.0 - gradyWeight[3];
                        gradyWeight[2] = 1.0 + gradyWeight[0] 
				- 2.0 * gradyWeight[3];
                        gradyWeight[1] = - gradyWeight[0] - gradyWeight[2] 
				- gradyWeight[3];
			break;
		case 4L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= (1.0 / 24.0) * xWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			xWeight[1] = t1 + t0;
			xWeight[3] = t1 - t0;
			xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
			xWeight[2] = 1.0 - xWeight[0] - xWeight[1] 
				- xWeight[3] - xWeight[4];
			/* grady */
                        w = y - (double)yIndex[2];
                        w2 = w * w;
                        gradyWeight[0] = - (1.0 / 8.0) * (w2 - w + 1.0 / 4.0);
                        t = (1.0 / 2.0) * w;
                        t0 = t * w - 11.0 / 24.0;
                        t1 = t - (2.0 / 3.0) * w2 * w;
                        gradyWeight[1] = t1 + t0;
                        gradyWeight[3] = t1 - t0;
                        gradyWeight[4] = gradyWeight[0] + t0 + (1.0 / 2.0);
                        gradyWeight[2] = - gradyWeight[0] - gradyWeight[1] 
				- gradyWeight[3] - gradyWeight[4];
			break;
		case 5L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- xWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			xWeight[2] = t0 + t1;
			xWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			xWeight[1] = t0 + t1;
			xWeight[4] = t0 - t1;
			/* grady */
                        w = y - (double)yIndex[2];
                        w2 = w * w;
                        gradyWeight[5] = (1.0 / 24.0) * w2 * w2;
                        w2 -= w;
                        w4 = w2 * w2;
                        w -= 1.0 / 2.0;
                        t = w2 * (w2 - 3.0);
                        gradyWeight[0] = (1.0 / 12.0) * w * (1.0 + 2.0 * w2) 
				- gradyWeight[5];
                        t0 = (1.0 / 12.0) * (2.0 * w2 - 5.0) * w;
                        t1 = (-1.0 / 12.0) * (w4 - 3.0 * w2 + 4.0) 
				- (1.0 / 6.0) * w * w * (2.0 * w2 - 3.0);
                        gradyWeight[2] = t0 + t1;
                        gradyWeight[3] = t0 - t1;
                        t0 = (1.0 / 8.0) * (3.0 - 2.0 * w2) * w;
                        t1 = (1.0 / 24.0) * (w4 - w2 - 5.0 
				+ 2.0 * w * w * (2.0 * w2 - 1));
                        gradyWeight[1] = t0 + t1;
                        gradyWeight[4] = t0 - t1;
			break;
		default:
			printf("Invalid spline degree\n");
			return(0.0);
	}

	/* apply the mirror boundary conditions */
	for (k = 0L; k <= SplineDegree; k++) {
		xIndex[k] = (Width == 1L) ? (0L) : ((xIndex[k] < 0L) ?
			(-xIndex[k] - Width2 * ((-xIndex[k]) / Width2))
			: (xIndex[k] - Width2 * (xIndex[k] / Width2)));
		if (Width <= xIndex[k]) {
			xIndex[k] = Width2 - xIndex[k];
		}
		yIndex[k] = (Height == 1L) ? (0L) : ((yIndex[k] < 0L) ?
			(-yIndex[k] - Height2 * ((-yIndex[k]) / Height2))
			: (yIndex[k] - Height2 * (yIndex[k] / Height2)));
		if (Height <= yIndex[k]) {
			yIndex[k] = Height2 - yIndex[k];
		}
	}

	/* perform interpolation */
	for (j = 0L; j <= SplineDegree; j++) {
		p = Bout + (ptrdiff_t)(yIndex[j] * Width);

		for (i = 0L; i <= SplineDegree; i++) {
			p[xIndex[i]]  += gradyWeight[j] * xWeight[i] * Bcoeff;
		}
	}

	return 1;

} /* end InterpolatedGradYTranMirr */

extern double	InterpolatedValueTranZero (
                        float   *Bout,      
                        long    Width,      
                        long    Height,    
                        float   Bcoeff,   
                        double  x,        
                        double  y,        
                        long    SplineDegree 
		)
{ /* begin InterpolatedValueTranZero */

	float	*p;
	double	xWeight[10], yWeight[10];
	double	w, w2, w4, t, t0, t1;
	long	xIndex[10], yIndex[10];
	long	i, j, k;

	/* compute the interpolation indexes */
	if (SplineDegree & 1L) {
		i = (long)floor(x) - SplineDegree / 2L;
		j = (long)floor(y) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
		}
	}
	else {
		i = (long)floor(x + 0.5) - SplineDegree / 2L;
		j = (long)floor(y + 0.5) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
		}
	}

	/* compute the interpolation weights */
	switch (SplineDegree) {
		case 2L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[1] = 3.0 / 4.0 - w * w;
			xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
			xWeight[0] = 1.0 - xWeight[1] - xWeight[2];
			/* y */
			w = y - (double)yIndex[1];
			yWeight[1] = 3.0 / 4.0 - w * w;
			yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
			yWeight[0] = 1.0 - yWeight[1] - yWeight[2];
			break;
		case 3L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[3] = (1.0 / 6.0) * w * w * w;
			xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- xWeight[3];
			xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
			xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];
			/* y */
			w = y - (double)yIndex[1];
			yWeight[3] = (1.0 / 6.0) * w * w * w;
			yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- yWeight[3];
			yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
			yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];
			break;
		case 4L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= (1.0 / 24.0) * xWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			xWeight[1] = t1 + t0;
			xWeight[3] = t1 - t0;
			xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
			xWeight[2] = 1.0 - xWeight[0] - xWeight[1] 
				- xWeight[3] - xWeight[4];
			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= (1.0 / 24.0) * yWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			yWeight[1] = t1 + t0;
			yWeight[3] = t1 - t0;
			yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
			yWeight[2] = 1.0 - yWeight[0] - yWeight[1] 
				- yWeight[3] - yWeight[4];
			break;
		case 5L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- xWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			xWeight[2] = t0 + t1;
			xWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			xWeight[1] = t0 + t1;
			xWeight[4] = t0 - t1;
			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- yWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			yWeight[2] = t0 + t1;
			yWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			yWeight[1] = t0 + t1;
			yWeight[4] = t0 - t1;
			break;
		case 6L:
			/* x */
			w = x - (double)xIndex[3];
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0] * xWeight[0];
			xWeight[0] *= xWeight[0] / 720.0;
			xWeight[1] = (361.0 / 192.0 - w * (59.0 / 8.0 + w
				* (-185.0 / 16.0 + w * (25.0 / 3.0 
				+ w * (-5.0 / 2.0 + w)
				* (1.0 / 2.0 + w))))) / 120.0;
			xWeight[2] = (10543.0 / 960.0 + w * (-289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (43.0 / 6.0 
				+ w * (-17.0 / 4.0 + w
				* (-1.0 + w)))))) / 48.0;
			w2 = w * w;
			xWeight[3] = (5887.0 / 320.0 - w2 * (231.0 / 16.0 - w2
				* (21.0 / 4.0 - w2))) / 36.0;
			xWeight[4] = (10543.0 / 960.0 + w * (289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (-43.0 / 6.0 
				+ w * (-17.0 / 4.0 + w
				* (1.0 + w)))))) / 48.0;
			xWeight[6] = 1.0 / 2.0 + w;
			xWeight[6] *= xWeight[6] * xWeight[6];
			xWeight[6] *= xWeight[6] / 720.0;
			xWeight[5] = 1.0 - xWeight[0] - xWeight[1] 
				- xWeight[2] - xWeight[3]
				- xWeight[4] - xWeight[6];
			/* y */
			w = y - (double)yIndex[3];
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0] * yWeight[0];
			yWeight[0] *= yWeight[0] / 720.0;
			yWeight[1] = (361.0 / 192.0 - w * (59.0 / 8.0 + w
				* (-185.0 / 16.0 + w * (25.0 / 3.0 
				+ w * (-5.0 / 2.0 + w)
				* (1.0 / 2.0 + w))))) / 120.0;
			yWeight[2] = (10543.0 / 960.0 + w * (-289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (43.0 / 6.0 
				+ w * (-17.0 / 4.0 + w
				* (-1.0 + w)))))) / 48.0;
			w2 = w * w;
			yWeight[3] = (5887.0 / 320.0 - w2 * (231.0 / 16.0 - w2
				* (21.0 / 4.0 - w2))) / 36.0;
			yWeight[4] = (10543.0 / 960.0 + w * (289.0 / 16.0 + w
				* (79.0 / 16.0 + w * (-43.0 / 6.0 
				+ w * (-17.0 / 4.0 + w * (1.0 + w)))))) / 48.0;
			yWeight[6] = 1.0 / 2.0 + w;
			yWeight[6] *= yWeight[6] * yWeight[6];
			yWeight[6] *= yWeight[6] / 720.0;
			yWeight[5] = 1.0 - yWeight[0] - yWeight[1] - yWeight[2] 
				- yWeight[3] - yWeight[4] - yWeight[6];
			break;
		case 7L:
			/* x */
			w = x - (double)xIndex[3];
			xWeight[0] = 1.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0] * xWeight[0];
			xWeight[0] *= (1.0 - w) / 5040.0;
			w2 = w * w;
			xWeight[1] = (120.0 / 7.0 + w * (-56.0 + w * (72.0 + w
				* (-40.0 + w2*(12.0 + w*(-6.0 + w)))))) / 720.0;
			xWeight[2] = (397.0/7.0 - w*(245.0/3.0 + w * (-15.0 + w
				* (-95.0 / 3.0 + w * (15.0 + w * (5.0 + w
				* (-5.0 + w))))))) / 240.0;
			xWeight[3] = (2416.0 / 35.0 + w2*(-48.0 + w2*(16.0 + w2
				* (-4.0 + w)))) / 144.0;
			xWeight[4] = (1191.0 / 35.0 - w * (-49.0 + w * (-9.0 + w
				* (19.0 + w*(-3.0 + w)*(-3.0 + w2))))) / 144.0;
			xWeight[5] = (40.0 / 7.0 + w*(56.0 / 3.0 + w * (24.0 + w
				*(40.0/3.0 + w2*(-4.0 + w*(-2.0 + w))))))/240.0;
			xWeight[7] = w2;
			xWeight[7] *= xWeight[7] * xWeight[7];
			xWeight[7] *= w / 5040.0;
			xWeight[6] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] 
				- xWeight[3]
				- xWeight[4] - xWeight[5] - xWeight[7];
			/* y */
			w = y - (double)yIndex[3];
			yWeight[0] = 1.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0] * yWeight[0];
			yWeight[0] *= (1.0 - w) / 5040.0;
			w2 = w * w;
			yWeight[1] = (120.0 / 7.0 + w * (-56.0 + w * (72.0 + w
				*(-40.0 + w2*(12.0 + w*(-6.0 + w)))))) / 720.0;
			yWeight[2] = (397.0 / 7.0 - w*(245.0/3.0 + w*(-15.0 + w
				* (-95.0 / 3.0 + w * (15.0 + w * (5.0 + w
				* (-5.0 + w))))))) / 240.0;
			yWeight[3] = (2416.0/35.0 + w2*(-48.0 + w2 * (16.0 + w2
				* (-4.0 + w)))) / 144.0;
			yWeight[4] = (1191.0 / 35.0 - w * (-49.0 + w * (-9.0 + w
				* (19.0 + w*(-3.0 + w)*(-3.0 + w2))))) / 144.0;
			yWeight[5] = (40.0/7.0 + w * (56.0 / 3.0 + w * (24.0 + w
				* (40.0 / 3.0 + w2 * (-4.0 + w * (-2.0 
				+ w)))))) / 240.0;
			yWeight[7] = w2;
			yWeight[7] *= yWeight[7] * yWeight[7];
			yWeight[7] *= w / 5040.0;
			yWeight[6] = 1.0 - yWeight[0] - yWeight[1] - yWeight[2] 
				- yWeight[3]
				- yWeight[4] - yWeight[5] - yWeight[7];
			break;
		case 8L:
			/* x */
			w = x - (double)xIndex[4];
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0] / 40320.0;
			w2 = w * w;
			xWeight[1] = (39.0/16.0 - w*(6.0 + w*(-9.0 / 2.0 + w2)))
				* (21.0 / 16.0 + w * (-15.0/4.0 + w*(9.0/2.0 + w
				* (-3.0 + w)))) / 5040.0;
			xWeight[2] = (82903.0 / 1792.0 + w * (-4177.0 / 32.0 + w
				* (2275.0 / 16.0 + w * (-487.0 / 8.0 
				+ w * (-85.0 / 8.0 + w
				* (41.0 / 2.0 + w * (-5.0 + w * (-2.0 
				+ w)))))))) / 1440.0;
			xWeight[3] = (310661.0/1792.0 - w * (14219.0 / 64.0 + w
				* (-199.0 / 8.0 + w * (-1327.0 / 16.0 
				+ w * (245.0 / 8.0 + w
				* (53.0 / 4.0 + w * (-8.0 + w * (-1.0 
				+ w)))))))) / 720.0;
			xWeight[4] = (2337507.0/8960.0 + w2*(-2601.0 / 16.0 + w2
				* (387.0 / 8.0 + w2 * (-9.0 + w2)))) / 576.0;
			xWeight[5] = (310661.0 / 1792.0 - w*(-14219.0 / 64.0 + w
				* (-199.0 / 8.0 + w * (1327.0 / 16.0 
				+ w * (245.0 / 8.0 + w
				* (-53.0 / 4.0 + w * (-8.0 + w * (1.0 
				+ w)))))))) / 720.0;
			xWeight[7] = (39.0/16.0 - w*(-6.0 + w*(-9.0/2.0 + w2)))
				* (21.0/16.0 + w*(15.0/4.0 + w * (9.0 / 2.0 + w
				* (3.0 + w)))) / 5040.0;
			xWeight[8] = 1.0 / 2.0 + w;
			xWeight[8] *= xWeight[8];
			xWeight[8] *= xWeight[8];
			xWeight[8] *= xWeight[8] / 40320.0;
			xWeight[6] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] 
				- xWeight[3] - xWeight[4] - xWeight[5] 
				- xWeight[7] - xWeight[8];
			/* y */
			w = y - (double)yIndex[4];
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0] / 40320.0;
			w2 = w * w;
			yWeight[1] = (39.0 / 16.0 - w*(6.0 + w*(-9.0/2.0 + w2)))
				* (21.0/16.0 + w*(-15.0/4.0 + w * (9.0 / 2.0 + w
				* (-3.0 + w)))) / 5040.0;
			yWeight[2] = (82903.0 / 1792.0 + w * (-4177.0 / 32.0 + w
				* (2275.0 / 16.0 + w * (-487.0 / 8.0 
				+ w * (-85.0 / 8.0 + w
				* (41.0 / 2.0 + w * (-5.0 + w * (-2.0 
				+ w)))))))) / 1440.0;
			yWeight[3] = (310661.0/1792.0 - w * (14219.0 / 64.0 + w
				* (-199.0 / 8.0 + w * (-1327.0 / 16.0 
				+ w * (245.0 / 8.0 + w
				* (53.0 / 4.0 + w * (-8.0 + w * (-1.0 
				+ w)))))))) / 720.0;
			yWeight[4] = (2337507.0/8960.0 + w2*(-2601.0 / 16.0 + w2
				* (387.0 / 8.0 + w2 * (-9.0 + w2)))) / 576.0;
			yWeight[5] = (310661.0 / 1792.0 - w*(-14219.0 / 64.0 + w
				* (-199.0 / 8.0 + w * (1327.0 / 16.0 
				+ w * (245.0 / 8.0 + w
				* (-53.0 / 4.0 + w * (-8.0 + w * (1.0 
				+ w)))))))) / 720.0;
			yWeight[7] = (39.0/16.0 - w*(-6.0 + w*(-9.0/2.0 + w2)))
				* (21.0 / 16.0 + w * (15.0/4.0 + w*(9.0/2.0 + w
				* (3.0 + w)))) / 5040.0;
			yWeight[8] = 1.0 / 2.0 + w;
			yWeight[8] *= yWeight[8];
			yWeight[8] *= yWeight[8];
			yWeight[8] *= yWeight[8] / 40320.0;
			yWeight[6] = 1.0 - yWeight[0] - yWeight[1] - yWeight[2] 
				- yWeight[3] - yWeight[4] - yWeight[5] 
				- yWeight[7] - yWeight[8];
			break;
		case 9L:
			/* x */
			w = x - (double)xIndex[4];
			xWeight[0] = 1.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0];
			xWeight[0] *= xWeight[0] * (1.0 - w) / 362880.0;
			xWeight[1] = (502.0 / 9.0 + w * (-246.0 + w * (472.0 + w
				* (-504.0 + w * (308.0 + w * (-84.0 
				+ w * (-56.0 / 3.0 + w
				* (24.0 + w * (-8.0 + w))))))))) / 40320.0;
			xWeight[2] = (3652.0 / 9.0 - w * (2023.0 / 2.0 
				+ w * (-952.0 + w
				* (938.0 / 3.0 + w * (112.0 + w * (-119.0 
				+ w * (56.0 / 3.0 + w
				* (14.0 + w * (-7.0 + w))))))))) / 10080.0;
			xWeight[3] = (44117.0 / 42.0 + w * (-2427.0 / 2.0 
				+ w * (66.0 + w
				* (434.0 + w * (-129.0 + w * (-69.0 
				+ w * (34.0 + w * (6.0 + w
				* (-6.0 + w))))))))) / 4320.0;
			w2 = w * w;
			xWeight[4] = (78095.0 / 63.0 - w2 * (700.0 
				+ w2 * (-190.0 + w2
				* (100.0 / 3.0 + w2 * (-5.0 + w))))) / 2880.0;
			xWeight[5] = (44117.0/63.0 + w * (809.0 + w * (44.0 + w
				* (-868.0 / 3.0 + w * (-86.0 + w * (46.0 
				+ w * (68.0 / 3.0 + w
				* (-4.0 + w * (-4.0 + w))))))))) / 2880.0;
			xWeight[6] = (3652.0 / 21.0 - w * (-867.0 / 2.0 
				+ w * (-408.0 + w
				* (-134.0 + w * (48.0 + w * (51.0 
				+ w * (-4.0 + w) * (-1.0 + w)
				* (2.0 + w))))))) / 4320.0;
			xWeight[7] = (251.0 / 18.0 + w * (123.0 / 2.0 
				+ w * (118.0 + w
				*(126.0 + w*(77.0 + w*(21.0 + w*(-14.0 / 3.0 + w
				* (-6.0 + w * (-2.0 + w))))))))) / 10080.0;
			xWeight[9] = w2 * w2;
			xWeight[9] *= xWeight[9] * w / 362880.0;
			xWeight[8] = 1.0 - xWeight[0] - xWeight[1] - xWeight[2] 
				- xWeight[3] - xWeight[4] - xWeight[5] 
				- xWeight[6] - xWeight[7] - xWeight[9];
			/* y */
			w = y - (double)yIndex[4];
			yWeight[0] = 1.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0];
			yWeight[0] *= yWeight[0] * (1.0 - w) / 362880.0;
			yWeight[1] = (502.0 / 9.0 + w * (-246.0 + w * (472.0 + w
				* (-504.0 + w * (308.0 + w * (-84.0 
				+ w * (-56.0 / 3.0 + w
				* (24.0 + w * (-8.0 + w))))))))) / 40320.0;
			yWeight[2] = (3652.0 / 9.0 - w * (2023.0 / 2.0 
				+ w * (-952.0 + w
				* (938.0 / 3.0 + w * (112.0 + w * (-119.0 
				+ w * (56.0 / 3.0 + w
				* (14.0 + w * (-7.0 + w))))))))) / 10080.0;
			yWeight[3] = (44117.0 / 42.0 + w * (-2427.0 / 2.0 
				+ w * (66.0 + w
				* (434.0 + w * (-129.0 + w * (-69.0 
				+ w * (34.0 + w * (6.0 + w
				* (-6.0 + w))))))))) / 4320.0;
			w2 = w * w;
			yWeight[4] = (78095.0/63.0 - w2*(700.0 + w2*(-190.0 + w2
				* (100.0 / 3.0 + w2 * (-5.0 + w))))) / 2880.0;
			yWeight[5] = (44117.0 / 63.0 + w * (809.0 + w*(44.0 + w
				* (-868.0 / 3.0 + w * (-86.0 + w * (46.0 
				+ w * (68.0 / 3.0 + w
				* (-4.0 + w * (-4.0 + w))))))))) / 2880.0;
			yWeight[6] = (3652.0 / 21.0 - w * (-867.0 / 2.0 
				+ w * (-408.0 + w
				* (-134.0 + w * (48.0 + w * (51.0 
				+ w * (-4.0 + w) * (-1.0 + w)
				* (2.0 + w))))))) / 4320.0;
			yWeight[7] = (251.0/18.0 + w*(123.0/2.0 + w * (118.0 + w
				* (126.0 + w*(77.0 + w*(21.0 + w*(-14.0/3.0 + w
				* (-6.0 + w * (-2.0 + w))))))))) / 10080.0;
			yWeight[9] = w2 * w2;
			yWeight[9] *= yWeight[9] * w / 362880.0;
			yWeight[8] = 1.0 - yWeight[0] - yWeight[1] - yWeight[2] 
				- yWeight[3] - yWeight[4] - yWeight[5] 
				- yWeight[6] - yWeight[7] - yWeight[9];
			break;
		default:
			printf("Invalid spline degree\n");
			return(0.0);
	}

	/* perform interpolation */
	for (j = 0L; j <= SplineDegree; j++) {
                if ((yIndex[j] >= 0) && (yIndex[j] < Height))
                {
                        p = Bout + (ptrdiff_t)(yIndex[j] * Width);

                        for (i = 0L; i <= SplineDegree; i++) {
                                if ((xIndex[i] >= 0) && (xIndex[i] < Width))
                                {
                                        p[xIndex[i]] += yWeight[j] * xWeight[i]
				 			* Bcoeff;
                                }
                        }
                }
        }

	return 1;

} /* end InterpolatedValueTranZero */

extern double	InterpolatedGradXTranZero (
                        float   *Bout,     
                        long    Width,   
                        long    Height,  
                        float   Bcoeff,     
                        double  x,       
                        double  y,       
                        long    SplineDegree  
		)
{ /* begin InterpolatedGradXTranZero */

	float	*p;
	double	gradxWeight[10], yWeight[10];
	double	w, w2, w4, t, t0, t1;
	long	xIndex[10], yIndex[10];
	long	Width2 = 2L * Width - 2L, Height2 = 2L * Height - 2L;
	long	i, j, k;

	/* compute the interpolation indexes */
	if (SplineDegree & 1L) {
		i = (long)floor(x) - SplineDegree / 2L;
		j = (long)floor(y) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
		}
	}
	else {
		i = (long)floor(x + 0.5) - SplineDegree / 2L;
		j = (long)floor(y + 0.5) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
		}
	}

	/* compute the interpolation weights */
	switch (SplineDegree) {
		case 2L:
			/* gradx */
                        w = x - (double)xIndex[1];
                        gradxWeight[1] = - 2.0 * w;
                        gradxWeight[2] = (1.0 / 2.0) * (1.0 - gradxWeight[1]);
                        gradxWeight[0] = - gradxWeight[1] - gradxWeight[2];
			/* y */
			w = y - (double)yIndex[1];
			yWeight[1] = 3.0 / 4.0 - w * w;
			yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
			yWeight[0] = 1.0 - yWeight[1] - yWeight[2];
			break;
		case 3L:
 			/* gradx */
                        w = x - (double)xIndex[1];
                        gradxWeight[3] = (1.0 / 2.0) * w * w;
                        gradxWeight[0] = w - 1.0 / 2.0 - gradxWeight[3];
                        gradxWeight[2] = 1.0 + gradxWeight[0] 
				- 2.0 * gradxWeight[3];
                        gradxWeight[1] = - gradxWeight[0] - gradxWeight[2] 					- gradxWeight[3];
			/* y */
			w = y - (double)yIndex[1];
			yWeight[3] = (1.0 / 6.0) * w * w * w;
			yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- yWeight[3];
			yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
			yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];
			break;
		case 4L:
 			/* gradx */
                        w = x - (double)xIndex[2];
                        w2 = w * w;
                        gradxWeight[0] = - (1.0 / 8.0) * (w2 - w + 1.0 / 4.0);
                        t = (1.0 / 2.0) * w;
                        t0 = t * w - 11.0 / 24.0;
                        t1 = t - (2.0 / 3.0) * w2 * w;
                        gradxWeight[1] = t1 + t0;
                        gradxWeight[3] = t1 - t0;
                        gradxWeight[4] = gradxWeight[0] + t0 + (1.0 / 2.0);
                        gradxWeight[2] = - gradxWeight[0] - gradxWeight[1] 
				- gradxWeight[3] - gradxWeight[4];
			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			yWeight[0] = 1.0 / 2.0 - w;
			yWeight[0] *= yWeight[0];
			yWeight[0] *= (1.0 / 24.0) * yWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			yWeight[1] = t1 + t0;
			yWeight[3] = t1 - t0;
			yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
			yWeight[2] = 1.0 - yWeight[0] - yWeight[1] 
				- yWeight[3] - yWeight[4];
			break;
		case 5L:
  			/* gradx */
                        w = x - (double)xIndex[2];
                        w2 = w * w;
                        gradxWeight[5] = (1.0 / 24.0) * w2 * w2;
                        w2 -= w;
                        w4 = w2 * w2;
                        w -= 1.0 / 2.0;
                        t = w2 * (w2 - 3.0);
                        gradxWeight[0] = (1.0 / 12.0) * w * (1.0 + 2.0 * w2) 
				- gradxWeight[5];
                        t0 = (1.0 / 12.0) * (2.0 * w2 - 5.0) * w;
                        t1 = (-1.0 / 12.0) * (w4 - 3.0 * w2 + 4.0) 
				- (1.0 / 6.0) * w * w * (2.0 * w2 - 3.0);
                        gradxWeight[2] = t0 + t1;
                        gradxWeight[3] = t0 - t1;
                        t0 = (1.0 / 8.0) * (3.0 - 2.0 * w2) * w;
                        t1 = (1.0 / 24.0) * (w4 - w2 - 5.0 
				+ 2.0 * w * w * (2.0 * w2 - 1));
                        gradxWeight[1] = t0 + t1;
                        gradxWeight[4] = t0 - t1;
			/* y */
			w = y - (double)yIndex[2];
			w2 = w * w;
			yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- yWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			yWeight[2] = t0 + t1;
			yWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			yWeight[1] = t0 + t1;
			yWeight[4] = t0 - t1;
			break;
		default:
			printf("Invalid spline degree\n");
			return(0.0);
	}

	/* perform interpolation */
	for (j = 0L; j <= SplineDegree; j++) {
                if ((yIndex[j] >= 0) && (yIndex[j] < Height))
                {
                        p = Bout + (ptrdiff_t)(yIndex[j] * Width);

                        for (i = 0L; i <= SplineDegree; i++) {
                                if ((xIndex[i] >= 0) && (xIndex[i] < Width))
                                {
                                        p[xIndex[i]] += yWeight[j] * 
							gradxWeight[i] * Bcoeff;
                                }
                        }
                }
        }

	return 1;

} /* end InterpolatedGradXTranZero */

extern double	InterpolatedGradYTranZero (
                        float   *Bout,     
                        long    Width,    
                        long    Height,   
                        float   Bcoeff,  
                        double  x,    
                        double  y,    
                        long    SplineDegree   
		)
{ /* begin InterpolatedGradYTranZero */

	float	*p;
	double	xWeight[10], gradyWeight[10];
	double	w, w2, w4, t, t0, t1;
	long	xIndex[10], yIndex[10];
	long	Width2 = 2L * Width - 2L, Height2 = 2L * Height - 2L;
	long	i, j, k;

	/* compute the interpolation indexes */
	if (SplineDegree & 1L) {
		i = (long)floor(x) - SplineDegree / 2L;
		j = (long)floor(y) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
		}
	}
	else {
		i = (long)floor(x + 0.5) - SplineDegree / 2L;
		j = (long)floor(y + 0.5) - SplineDegree / 2L;
		for (k = 0L; k <= SplineDegree; k++) {
			xIndex[k] = i++;
			yIndex[k] = j++;
		}
	}

	/* compute the interpolation weights */
	switch (SplineDegree) {
		case 2L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[1] = 3.0 / 4.0 - w * w;
			xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
			xWeight[0] = 1.0 - xWeight[1] - xWeight[2];
 			/* grady */
                        w = y - (double)yIndex[1];
                        gradyWeight[1] = - 2.0 * w;
                        gradyWeight[2] = (1.0 / 2.0) * (1.0 - gradyWeight[1]);
                        gradyWeight[0] = - gradyWeight[1] - gradyWeight[2];
			break;
		case 3L:
			/* x */
			w = x - (double)xIndex[1];
			xWeight[3] = (1.0 / 6.0) * w * w * w;
			xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) 
				- xWeight[3];
			xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
			xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];
			/* grady */
                        w = y - (double)yIndex[1];
                        gradyWeight[3] = (1.0 / 2.0) * w * w;
                        gradyWeight[0] = w - 1.0 / 2.0 - gradyWeight[3];
                        gradyWeight[2] = 1.0 + gradyWeight[0] 
				- 2.0 * gradyWeight[3];
                        gradyWeight[1] = - gradyWeight[0] - gradyWeight[2] 
				- gradyWeight[3];
			break;
		case 4L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			t = (1.0 / 6.0) * w2;
			xWeight[0] = 1.0 / 2.0 - w;
			xWeight[0] *= xWeight[0];
			xWeight[0] *= (1.0 / 24.0) * xWeight[0];
			t0 = w * (t - 11.0 / 24.0);
			t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
			xWeight[1] = t1 + t0;
			xWeight[3] = t1 - t0;
			xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
			xWeight[2] = 1.0 - xWeight[0] - xWeight[1] 
				- xWeight[3] - xWeight[4];
			/* grady */
                        w = y - (double)yIndex[2];
                        w2 = w * w;
                        gradyWeight[0] = - (1.0 / 8.0) * (w2 - w + 1.0 / 4.0);
                        t = (1.0 / 2.0) * w;
                        t0 = t * w - 11.0 / 24.0;
                        t1 = t - (2.0 / 3.0) * w2 * w;
                        gradyWeight[1] = t1 + t0;
                        gradyWeight[3] = t1 - t0;
                        gradyWeight[4] = gradyWeight[0] + t0 + (1.0 / 2.0);
                        gradyWeight[2] = - gradyWeight[0] - gradyWeight[1] 
				- gradyWeight[3] - gradyWeight[4];
			break;
		case 5L:
			/* x */
			w = x - (double)xIndex[2];
			w2 = w * w;
			xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
			w2 -= w;
			w4 = w2 * w2;
			w -= 1.0 / 2.0;
			t = w2 * (w2 - 3.0);
			xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) 
				- xWeight[5];
			t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
			t1 = (-1.0 / 12.0) * w * (t + 4.0);
			xWeight[2] = t0 + t1;
			xWeight[3] = t0 - t1;
			t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
			t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
			xWeight[1] = t0 + t1;
			xWeight[4] = t0 - t1;
			/* grady */
                        w = y - (double)yIndex[2];
                        w2 = w * w;
                        gradyWeight[5] = (1.0 / 24.0) * w2 * w2;
                        w2 -= w;
                        w4 = w2 * w2;
                        w -= 1.0 / 2.0;
                        t = w2 * (w2 - 3.0);
                        gradyWeight[0] = (1.0 / 12.0) * w * (1.0 + 2.0 * w2) 
				- gradyWeight[5];
                        t0 = (1.0 / 12.0) * (2.0 * w2 - 5.0) * w;
                        t1 = (-1.0 / 12.0) * (w4 - 3.0 * w2 + 4.0) 
				- (1.0 / 6.0) * w * w * (2.0 * w2 - 3.0);
                        gradyWeight[2] = t0 + t1;
                        gradyWeight[3] = t0 - t1;
                        t0 = (1.0 / 8.0) * (3.0 - 2.0 * w2) * w;
                        t1 = (1.0 / 24.0) * (w4 - w2 - 5.0 
				+ 2.0 * w * w * (2.0 * w2 - 1));
                        gradyWeight[1] = t0 + t1;
                        gradyWeight[4] = t0 - t1;
			break;
		default:
			printf("Invalid spline degree\n");
			return(0.0);
	}

	/* perform interpolation */
	for (j = 0L; j <= SplineDegree; j++) {
                if ((yIndex[j] >= 0) && (yIndex[j] < Height))
                {
                        p = Bout + (ptrdiff_t)(yIndex[j] * Width);

                        for (i = 0L; i <= SplineDegree; i++) {
                                if ((xIndex[i] >= 0) && (xIndex[i] < Width))
                                {
                                        p[xIndex[i]] += gradyWeight[j] * 
							xWeight[i] * Bcoeff;
                                }
                        }
                }
        }

	return 1;

} /* end InterpolatedGradYTranZero */
