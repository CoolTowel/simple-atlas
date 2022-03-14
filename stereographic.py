import numpy as np
from numpy import ma
from matplotlib import scale as mscale
from matplotlib import transforms as mtransforms
from matplotlib.ticker import FixedLocator, FuncFormatter
import matplotlib.font_manager as fm

font_path = './fonts/Libertinus-7.040/static/OTF/'
font_files = fm.findSystemFonts(fontpaths=font_path)
for font_file in font_files:
    fm.fontManager.addfont(font_file)
        
class StereographicZenithScale(mscale.ScaleBase):
    """
    Scales data in range -pi/2 to pi/2 (-90 to 90 degrees) using
    the system used to scale latitudes in a Mercator__ projection.

    The scale function:
      tan((pi/2-theta)/2)

    The inverse scale function:
      pi/2-2*(arctan(theta))

    Only show the sky above the horizon
    """

    # The scale class must have a member ``name`` that defines the string used
    # to select the scale.  For example, ``ax.set_yscale("mercator")`` would be
    # used to select this scale.
    name = 'stereographiczenith' # prevent confilt with astropy 

    def __init__(self, axis, *, thresh=0, **kwargs):
        """
        Any keyword arguments passed to ``set_xscale`` and ``set_yscale`` will
        be passed along to the scale's constructor.

        thresh: The degree less which to crop the data.
        """
        super().__init__(axis)
        if thresh >= np.pi / 2:
            raise ValueError("thresh must be less than pi/2")
        self.thresh = thresh

    def get_transform(self):
        """
        Override this method to return a new instance that does the
        actual transformation of the data.

        The StereographicLatitudeScale class is defined below as a
        nested class of this one.
        """
        return self.StereographicZenithTransform(self.thresh)

    def set_default_locators_and_formatters(self, axis):
        """
        Override to set up the locators and formatters to use with the
        scale.  This is only required if the scale requires custom
        locators and formatters.  Writing custom locators and
        formatters is rather outside the scope of this example, but
        there are many helpful examples in :mod:`.ticker`.

        In our case, the Mercator example uses a fixed locator from -90 to 90
        degrees and a custom formatter to convert the radians to degrees and
        put a degree symbol after the value.
        """
        fmt = FuncFormatter(
            lambda x, pos=None: f"{(90-np.degrees(x)):.0f}\N{DEGREE SIGN}")
        axis.set(major_locator=FixedLocator(np.radians([30,60])),
                 major_formatter=fmt, minor_formatter=fmt)

    def limit_range_for_scale(self, vmin, vmax, minpos):
        """
        Override to limit the bounds of the axis to the domain of the
        transform.  In the case of Mercator, the bounds should be
        limited to the threshold that was passed in.  Unlike the
        autoscaling provided by the tick locators, this range limiting
        will always be adhered to, whether the axis range is set
        manually, determined automatically or changed through panning
        and zooming.
        """
        return 0,np.pi/2 #max(vmin, -self.thresh), min(vmax, self.thresh)

    class StereographicZenithTransform(mtransforms.Transform):
        # There are two value members that must be defined.
        # ``input_dims`` and ``output_dims`` specify number of input
        # dimensions and output dimensions to the transformation.
        # These are used by the transformation framework to do some
        # error checking and prevent incompatible transformations from
        # being connected together.  When defining transforms for a
        # scale, which are, by definition, separable and have only one
        # dimension, these members should always be set to 1.
        input_dims = output_dims = 1

        def __init__(self, thresh):
            mtransforms.Transform.__init__(self)
            self.thresh = thresh

        def transform_non_affine(self, a):
            """
            This transform takes a numpy array and returns a transformed copy.
            Since the range of the Mercator scale is limited by the
            user-specified threshold, the input array must be masked to
            contain only valid values.  Matplotlib will handle masked arrays
            and remove the out-of-range data from the plot.  However, the
            returned array *must* have the same shape as the input array, since
            these values need to remain synchronized with values in the other
            dimension.
            """
            masked = ma.masked_where((a < self.thresh), a)
            if masked.mask.any():
                return ma.tan(masked/2)
            else:
                return np.tan(a/2)
        def inverted(self):
            """
            Override this method so Matplotlib knows how to get the
            inverse transform for this transform.
            """
            return StereographicZenithScale.InvertedStereographicZenithTransform(
                self.thresh)

    class InvertedStereographicZenithTransform(mtransforms.Transform):
        input_dims = output_dims = 1

        def __init__(self, thresh):
            mtransforms.Transform.__init__(self)
            self.thresh = thresh

        def transform_non_affine(self, a):
            return 2*np.arctan(a)

        def inverted(self):
            return StereographicZenithScale.InvertedStereographicZenithTransform(self.thresh)

# class StereographicAzimuthScale(mscale.ScaleBase):
#     """
#     The scale function:
#       theta+90

#     The inverse scale function:
#       theta-90
#     """

#     # The scale class must have a member ``name`` that defines the string used
#     # to select the scale.  For example, ``ax.set_yscale("mercator")`` would be
#     # used to select this scale.
#     name = 'azimuth' # prevent confilt with astropy 

#     def __init__(self, axis, **kwargs):
#         """
#         Any keyword arguments passed to ``set_xscale`` and ``set_yscale`` will
#         be passed along to the scale's constructor.
#         """
#         super().__init__(axis)

#     def get_transform(self):
#         """
#         Override this method to return a new instance that does the
#         actual transformation of the data.
#         The StereographicAzimuthScale class is defined below as a
#         nested class of this one.
#         """
#         return self.StereographicAzimuthTransform()
        
#     def set_default_locators_and_formatters(self, axis):
#         """
#         Override to set up the locators and formatters to use with the
#         scale.  This is only required if the scale requires custom
#         locators and formatters.  Writing custom locators and
#         formatters is rather outside the scope of this example, but
#         there are many helpful examples in :mod:`.ticker`.
#         """
#         fmt = FuncFormatter(azimuth_name)
#         axis.set(major_locator=FixedLocator([0,np.pi/2,np.pi,np.pi/2*3]),
#                  major_formatter=fmt, minor_formatter=fmt)

#     # def limit_range_for_scale(self, vmin, vmax, minpos):
#     #     """
#     #     Override to limit the bounds of the axis to the domain of the
#     #     transform.  In the case of Mercator, the bounds should be
#     #     limited to the threshold that was passed in.  Unlike the
#     #     autoscaling provided by the tick locators, this range limiting
#     #     will always be adhered to, whether the axis range is set
#     #     manually, determined automatically or changed through panning
#     #     and zooming.
#     #     """
#     #     return 0,np.pi/2 #max(vmin, -self.thresh), min(vmax, self.thresh)

#     class StereographicAzimuthTransform(mtransforms.Transform):
#         # There are two value members that must be defined.
#         # ``input_dims`` and ``output_dims`` specify number of input
#         # dimensions and output dimensions to the transformation.
#         # These are used by the transformation framework to do some
#         # error checking and prevent incompatible transformations from
#         # being connected together.  When defining transforms for a
#         # scale, which are, by definition, separable and have only one
#         # dimension, these members should always be set to 1.
#         input_dims = output_dims = 1

#         def __init__(self):
#             mtransforms.Transform.__init__(self)

#         def transform_non_affine(self, a):
#             """
#             This transform takes a numpy array and returns a transformed copy.
#             Since the range of the Mercator scale is limited by the
#             user-specified threshold, the input array must be masked to
#             contain only valid values.  Matplotlib will handle masked arrays
#             and remove the out-of-range data from the plot.  However, the
#             returned array *must* have the same shape as the input array, since
#             these values need to remain synchronized with values in the other
#             dimension.
#             """
#             return a
        
#         def inverted(self):
#             """
#             Override this method so Matplotlib knows how to get the
#             inverse transform for this transform.
#             """
#             return StereographicAzimuthScale.InvertedStereographicLatitudeTransform()

#     class InvertedStereographicLatitudeTransform(mtransforms.Transform):
#         input_dims = output_dims = 1

#         def __init__(self):
#             mtransforms.Transform.__init__(self)

#         def transform_non_affine(self, a):
#             return a

#         def inverted(self):
#             return StereographicAzimuthScale.InvertedStereographicLatitudeTransform()
        

# Now that the Scale class has been defined, it must be registered so
# that Matplotlib can find it.
mscale.register_scale(StereographicZenithScale)
# mscale.register_scale(StereographicAzimuthScale)



if __name__ == '__main__':
    import matplotlib.pyplot as plt

    r = np.linspace(30, 60, 50)/180*np.pi
    theta = np.linspace(0, 2 * np.pi, 50)

    plt.figure()
    plt.subplot(111, projection="polar")
    plt.plot(theta,r, '-', lw=2,)
    ax = plt.gca()
    ax.set_rscale('stereographiczenith')
    plt.xscale('azimuth')
    
    ax = plt.gca()
    ax.set_theta_zero_location("N")
    plt.grid(True)

    plt.show()