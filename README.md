# Stress-Testing Chatterjee's Xi Correlation Coefficient

Paquete de artículo + reproducibilidad para un estudio de simulación
Monte Carlo que somete a prueba el coeficiente de correlación xi de
Chatterjee (Chatterjee, 2021) contra cinco factores fundamentados
teóricamente: tamaño de muestra, tipo de variable, forma de la
distribución, estructura de dependencia, y régimen de dependencia.
Manuscrito preparado para envío a la Austrian Journal of Statistics (AJS).

## Estructura del repositorio

- `paper/` — fuente LaTeX del manuscrito (estilo AJS/JSS), bibliografía,
  y clase/estilo oficiales de la revista.
- `simulation/` — scripts de R para la simulación (mecanismos
  generadores de datos, calibración del régimen de dependencia,
  cómputo de los estadísticos, generación de tablas), resultados, y
  README propio con instrucciones de ejecución y semilla fija — sigue
  el mismo patrón de paquete reproducible usado en las validaciones
  Monte Carlo de AssumptionsLab.

## Estado

Trabajo en curso — el manuscrito y la simulación se están desarrollando
en paralelo. Ver `simulation/README.md` para el estado específico del
código de simulación.
