Fast rotate of images

fastrotate_yaro.m	Rotation code as provided by Yaroslavsky. The code works only for even-sized images, and requires the functions fftnorm and ifftnorm. The function seems to fail for angles that are around 180.

fastrotate_ref.m	Adaptation of Yaroslavsky's fastrotate_yaro to work with both even and odd sizes. The functions also does not require the functions fftnorm and ifftnorm. Since does not work either for angles around 180 degrees.

fastrotate.m		Optimized implementation of fastrotate_ref.m
