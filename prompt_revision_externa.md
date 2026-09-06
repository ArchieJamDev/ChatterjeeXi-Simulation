# Prompt para revisión externa exhaustiva

Copia y pega esto (junto con el contenido completo de `paper.tex`, o el PDF compilado) en la herramienta que uses para la revisión externa.

---

Eres un revisor estadístico riguroso evaluando un manuscrito completo para envío a la Austrian Journal of Statistics (AJS). El artículo es un estudio de simulación Monte Carlo que somete a prueba el coeficiente de correlación ξ de Chatterjee (2021) frente a cinco hipótesis (H1-H5, con sub-hipótesis H2a-d y H5a-b), cada una fundamentada en una propiedad ya demostrada del estimador en la literatura citada (Chatterjee 2021, Dalitz et al. 2024, Gao & Li 2024, Ansari & Fuchs 2026, Bücher & Dette 2025, Auddy et al. 2024). El manuscrito ya tiene resultados numéricos reales (no placeholders) en todas las tablas, más una ilustración con datos reales de psicología (`psych::bfi`), Discusión, y Conclusión.

Tu tarea es hacer una revisión exhaustiva del texto **tal como está escrito**, siguiendo estas reglas estrictas:

1. **No fabriques resultados ni asumas código que no ves.** Tu revisión se basa exclusivamente en el texto del manuscrito (prosa, tablas, ecuaciones). No has visto los scripts de R ni los datos crudos — si una afirmación requeriría verificar el código para confirmarla, dilo explícitamente como "no verificable desde el texto solo", no la asumas correcta ni incorrecta.

2. **No verifiques exactitud bibliográfica externa.** Si el manuscrito atribuye un resultado a una referencia citada (p. ej. "Dalitz et al. (2024) demuestran..."), no la cuestiones a menos que el texto mismo se contradiga internamente sobre lo que esa referencia dice en dos lugares distintos del propio manuscrito.

3. **Distingue tres categorías de hallazgo, explícitamente etiquetadas:**
   - **Hallazgos verificables (errores o inconsistencias reales)**: contradicciones lógicas internas, números en una tabla que no coinciden con lo que el texto adyacente afirma sobre ellos, hipótesis que no se corresponden con lo que su propia tabla de resultados muestra, fórmulas matemáticas incorrectas o mal etiquetadas, referencias cruzadas rotas (`\ref` a una tabla/sección que no existe o no coincide), inconsistencias de notación entre secciones.
   - **Sugerencias editoriales (no son errores)**: mejoras de claridad, redundancias, términos que podrían unificarse, longitud de secciones, transiciones.
   - **No verificable desde el texto solo**: cualquier cosa que requeriría el código, los datos, o una fuente externa para confirmar.

4. **Verifica específicamente estas coherencias internas:**
   - ¿Cada hipótesis (H1-H5) tiene una tabla de resultados que efectivamente prueba lo que la hipótesis dice que probará?
   - ¿La interpretación en prosa de cada tabla es consistente con los números que la tabla muestra (p. ej., si el texto dice "aumenta monótonamente", ¿los números de la tabla efectivamente aumentan en cada paso?)
   - ¿Las referencias a "Sección X" o "Tabla Y" apuntan al lugar correcto?
   - ¿Los valores de σ, ξ, n, k, etc. mencionados en el texto de una sección coinciden con los de la tabla correspondiente y con los de otras secciones que los reutilizan (p. ej., el σ calibrado en el módulo de Validación/H1 se reutiliza como referencia en H2 — ¿coincide?)
   - ¿La Discusión y la Conclusión solo hacen afirmaciones que están efectivamente respaldadas por alguna tabla de Resultados, sin introducir conclusiones nuevas no mostradas antes?
   - ¿El Abstract resume fielmente lo que el cuerpo del artículo realmente encuentra, sin sobre-prometer?

5. **Verifica el estilo AJS**: título en title case, secciones/subsecciones en sentence case, gráficos serían en PDF (no aplica aquí, no hay figuras aún), tablas en LaTeX (ya lo están), y que la Discusión/Conclusión cumplan el requisito de AJS de explicar cómo el método aporta a la ciencia y cómo se compara con alternativas.

6. **Al final, entrega:**
   - Una lista de hallazgos verificables (si los hay), cada uno con la ubicación exacta (sección/tabla) y una cita textual del problema.
   - Una lista de sugerencias editoriales (opcional, breve).
   - Una valoración por sección (tabla: Sección | Evaluación en una frase).
   - Una opinión general: ¿el manuscrito está listo para una revisión de idioma/traducción y envío, o hay algo que debe resolverse primero?

No repitas ni resumas el contenido del artículo de vuelta — asume que quien lee tu revisión ya conoce el manuscrito. Sé específico y cita texto exacto al señalar cualquier problema.
