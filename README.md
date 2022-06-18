# Mitsuba2 Plenoptic Light-Transport Probing
This repogitory contains the extensions of [Mitsuba2 renderer](https://github.com/mitsuba-renderer/mitsuba2) for light-transport rendering in consideration of light dimensions: polarization, direction, and time. 
We developed this software for our SIGGRAPH Asia 2021 publication: [Polarimetric Spatio-Temporal Light Transport Probing](https://light.princeton.edu/publication/polarimetric-spatio-temporal-light-transport-probing/). 

## Goal
We aim to simulate 7D light transport $\mathcal{T}$: two dimensions for illumination's spatial axes, two dimensions for camera's spatial axes, one for temporal dimension, and two for polarimetric changes. 

## Spatial dimensions
To simulate spatial imaging, we place a 2D emitter and a 2D detector and pixels are sampled in their screen spaces.
That is, we sample $uv\in[0,1]^2\sim\mathbb{S}_{\text{sensor}}, st\in[0,1]^2\sim\mathbb{S}_{\text{emitter}}$ as emitter and detector coordinates, and we only consider those that satisfy $\delta-\epsilon < st - uv < \delta$, where $\delta\in[0,1]^2,\epsilon\in[0,2]^2$. This means that we can specify a certain range of pixels on the light source to illuminate the scene, and we can specify specific pixel on the camera to sample the transport.

Use the following xml commands for setting the parameters of the spatial probing.

```
  <!-- Some measure of spatial probing on direct light contribution. For example, if we sample a
  light towards a given light towards the current point, it should only contribute if its own
  uv coordinate is within some range of uv. -->
  <float name="delta_x" value="$delta_x"/>
  <float name="eps_x" value="$eps_x"/>

  <!-- isolate light emitted between two y bounds given in 0 to 1. -->
  <!-- as compared to delta_x or delta_y and eps_x or eps_y above, this is directly on the
  emitter -->
  <float name="min_isolate_y" value="$min_isolate_y"/>
  <float name="max_isolate_y" value="$max_isolate_y"/>
```

## Temporal dimension
Temporal rendering is based on the time-of-flight from the emitter to the sensor. Only photons that satisfy $t_{\text{min}} < t_{\text{time of flight}} < t_{\text{max}}$ are rendered. This gating constraint can be used to only illuminate surfaces at certain summed distance between the camera and the emitter.

Use the following xml commands for setting the parameters of the temporal probing.

```
  <!-- time gated values, along with the speed to use when computing distance to time. -->
  <float name="min_t" value="$min_t"/>
  <float name="max_t" value="$max_t"/>
  <float name="light_speed" value="$light_speed"/>

  <!-- What time value will actually get rendered i.e. below is at 0.1 units. -->
  <float name="render_at" value="0.1"/>
<!-- The speed of light relative to distance. i.e. 1 distance unit is 10k time units. -->
  <float name="speed_of_light" value="10000."/>

```


## Polarimetric dimension
We make use of polarimetric-rendering feature implemented in Mitsuba2. Specifically, we constructed a virtual imaging system with dual-rotating retarders (see our paper) in Mitsuba2. 


## Usage
You can utilize our implemented features by simply writing your Mitsuba2 xml files.

```
<!-- A modified volumetric path integrator, which also allows for gating. -->
<integrator type="volpath_simple">
  <!-- time gated values, along with the speed to use when computing distance to time. -->
  <float name="min_t" value="$min_t"/>
  <float name="max_t" value="$max_t"/>
  <float name="light_speed" value="$light_speed"/>

  <!-- Some measure of spatial probing on direct light contribution. For example, if we sample a
  light towards a given light towards the current point, it should only contribute if its own
  uv coordinate is within some range of uv. -->
  <float name="delta_x" value="$delta_x"/>
  <float name="eps_x" value="$eps_x"/>

  <!-- whether to render things within give time, or render everything outside of given time -->
  <boolean name="invert_x" value="$invert_x"/>

  <!-- max number of bounces -->
  <integer name="max_depth" value="$max_depth"/>

  <!-- isolate light emitted between two y bounds given in 0 to 1-->
  <!-- as compared to delta_x or delta_y and eps_x or eps_y above, this is directly on the
  emitter -->
  <float name="min_isolate_y" value="$min_isolate_y"/>
  <float name="max_isolate_y" value="$max_isolate_y"/>
</integrator>
```

## About
Our code is based on the original [Mitsuba2 renderer](https://github.com/mitsuba-renderer/mitsuba2). [Julian Knodt](https://julianknodt.github.io/) mainly developed the extensions and [Seung-Hwan Baek](https://www.shbaek.com/) advised the process. If you use our renderer in your research, please cite our paper using the following BibText template:


```bib
@Article{baek2021polarLT,
   Author = {Seung-Hwan Baek and Felix Heide},
   Year = {2021},
   Title = {Polarimetric Spatio-Temporal Light Transport Probing},
   journal = {ACM Transactions on Graphics (Prof. SIGGRAPH Asia 2021)},
}
```