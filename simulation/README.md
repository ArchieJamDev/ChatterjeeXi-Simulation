# Paquete de reproducibilidad — simulación

Pendiente de implementación. Esta carpeta contendrá, siguiendo el mismo
patrón usado en las validaciones Monte Carlo de AssumptionsLab (script +
resultados + README por corrida, con semilla fija):

- Scripts de generación de datos para las cinco estructuras de
  dependencia (Sección 3.2 del manuscrito: lineal, cuadrática,
  senoidal, circular, heterocedástica).
- Script de calibración numérica del nivel de ruido `sigma` para
  alcanzar cada régimen objetivo de xi poblacional (H5).
- Script de discretización para las variantes binaria/ordinal (H2).
- Cómputo de xi_n, dCor, y D de Hoeffding / tau* sobre cada condición
  simulada.
- Resultados crudos (CSV) y scripts que generan las Tablas 2-6 del
  manuscrito directamente a partir de esos resultados.

## Software

Pendiente de especificar (versión de R, paquetes: `XICOR`, `energy`,
etc. — ver Sección 3.4 del manuscrito).

## Semilla

Pendiente de fijar antes de la primera corrida real.
