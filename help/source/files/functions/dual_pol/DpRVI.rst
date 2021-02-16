DpRVI
===================
Dual-pol Radar Vegetation Index .

.. math::

    \text{DpRVI} = 1- \text{DOP}_{dp}\Big(\frac{\lambda_1}{\lambda_1+\lambda_2}\Big)

where,

.. math::
	\text{DOP}_{dp} = \sqrt{1-\frac{4\times \text{det ([C2])}}{\text{(Trace [C2])}^2}}

:math:`\text{[C2]}` is co-variance matrix,  and :math:`\lambda_1, \lambda_2` are the eigen values of :math:`\langle\mathbf{[C2]}\rangle` matrix in descending order.

