Introduction
=============

KiMoPack is a project for the handling of spectral data measure at
multiple time-points. The current design is optimised for the use with
optical transient absorption data, but it has been successfully adapted
for the use with transient x-ray emission and spectro-electro chemistry
data.

It focuses on the main tasks an experimentator has
Loading and shaping of experiments, plotting of experiments, comparing of experiments,
analysing experiments with fast and/or advanced fitting routines and saving/exporting/presenting 
the results. 

The software can be used on several different levels. The simplest level packs everything 
into an object "TA" that contains all the parameters that are typically set. 
These objects also contain the typical functions that are used in an analysis. 
See :ref:`Main Tasks overview` for an overview of these functions. 
All active functions have a capital letter in the beginning.

At the lower levels a series of convenience functions for the efficient plotting of
one or two dimensional data is provided. These are typical in the main module 

For typical use a series of juypter notebooks are provided that guide 
through the a number of different use scenarios, and are suggesting the 
parameter that are typically set.

In addition a series of tutorial notebooks are provided that guide the user through the different functions. These Tutorials can either be downloaded or executed on a "mybinder" server via this badge.
 .. image:: https://mybinder.org/badge_logo.svg		  
	:target: https://mybinder.org/v2/gh/erdzeichen/KiMoPack/HEAD
	
In addition a small series of videos were produced to introduce the features and usage of KiMoPack: https://www.youtube.com/channel/UCmhiK0P9wXXjs_PJaitx8BQ
