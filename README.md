# Lab 14
## Relativistic particle motion


In this lab we will consider the motion of a charged particle which is under the influence of a strong pulsed electromagnetic field. The field is assumed to vary spatially only along the *x* direction. Furthermore, we assume that the only component of the electric field is the *y* component *Ey*. The only component of the magnetic field is then *Bz*.

To integrate the equation of motion

<p align="center">
<img src="stuffy_stuff/f1.png" width="200">
</p>

we use the **Boris push**, which is a splitting method that splits the motion due to the electric field from that due to the magnetic field. The update of the particle position is

<p align="center">
<img src="stuffy_stuff/f2.png" width="200">
</p>

The update of the particle momentum is

<p align="center">
<img src="stuffy_stuff/f3.png" width="550">
</p>

---

The source code contains a class `particle` which contains interfaces
to several member functions that you have to implement. Even though for our test problem the electric and the magnetic field have only a single component, implement the full three-dimensional push.
