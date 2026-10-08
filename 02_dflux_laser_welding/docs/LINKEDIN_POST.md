# Proposta de Publicação para o LinkedIn / LinkedIn Post Proposal

Aqui encontras duas versões prontas a publicar: uma em **Português** e outra em **Inglês**. Podes ajustar os links do teu GitHub ou adicionar uma imagem/vídeo dos teus resultados de simulação no Abaqus Viewer.

---

## Opção A: Português (Foco em Engenharia Térmica, Metalurgia & Simulação FEA)

🔥 **Modelagem Térmica de Soldadura Laser Circular em Abaqus: Subrotinas Fortran, Cinética de Fases e Avaliação de BTR**

A simulação de processos de soldadura a laser (laser beam welding) em componentes mecânicos de alta responsabilidade exige mais do que apenas aplicar uma condição de contorno térmica convencional. O balanço de energia, a cinemática da fonte e a evolução microestrutural na zona termicamente afetada (ZTA) são cruciais para prever tensões residuais e o risco de fissuração a quente.

Partilho um projeto e repositório de referência que desenvolvi para este efeito:

🔹 **Subrotina DFLUX em Fortran (F77 estrito e thread-safe):**
Implementação de uma fonte volumétrica cónica com distribuição gaussiana radial adaptativa $q(r, z) = Q_0 \exp(-3 r^2 / R_0(z)^2)$ em movimento circular contínuo, com controlo angular de rampas de subida, patamar e descida (para evitar defeitos no fecho do cordão). A normalização analítica garante 100% de conservação da potência absorvida $\eta Q_{\text{tot}}$.

🔹 **Metalurgia Avançada com `ABQ_PHASE_TRANS`:**
Integração do módulo metalúrgico do Abaqus/Standard para rastreio em tempo real de:
- Fusão e calor latente de solidificação.
- Cinética difusional (Austenite $\to$ Ferrite/Bainite) via formulação de Johnson-Mehl-Avrami (JMA).
- Transformação displaciva até Martensite via Koistinen-Marburger.
- Estado RLS (Raw-Liquid-Solid) calibrado para materiais sólidos forjados/usinados.

🔹 **Estratégia de Passo Temporal para BTR (Brittleness Temperature Range):**
Na soldadura a laser, o intervalo de fragilidade térmica (BTR, entre liquidus e solidus) é o momento crítico para nucleação de fissuras de solidificação. Para capturar com rigor esta transição, o primeiro passo de cálculo mantém um incremento de tempo estritamente controlado e uniforme durante e imediatamente após a extinção do feixe, transitando apenas no passo de arrefecimento para stepping adaptativo (`DELTMX = 25 °C`).

🔹 **Verificação Pré-Solver em Python:**
Antes de enviar o cálculo para o solver, uma rotina Python com integração numérica de Gauss $2\times 2\times 2$ nos elementos DC3D8 valida o balanço energético da malha discretizada, garantindo menos de 1.9% de erro de discretização.

Código e documentação completa disponíveis no GitHub: [link para o teu repositório]

#Abaqus #FEA #Welding #LaserWelding #Fortran #Python #MaterialsEngineering #Simulation #MechanicalEngineering #Metallurgy

---

## Opção B: English (Engineering & Computational Mechanics Focus)

🔥 **Simulating Circumferential Laser Welding in Abaqus: Fortran DFLUX, Phase Transformations & BTR Tracking**

Predicting residual stresses, distortion, and hot cracking susceptibility in laser-welded axisymmetric powertrain assemblies requires rigorous coupling between beam kinematics, energy conservation, and metallurgical phase transformations.

I have assembled an open benchmark repository demonstrating an end-to-end simulation methodology in Abaqus/Standard:

Key Highlights:
1. **Thread-Safe Fortran DFLUX Subroutine:**
A moving conical Gaussian volumetric heat source along a circular trajectory with dedicated power ramp schedules. The peak flux is analytically normalized over the conical geometry, ensuring exact theoretical energy input ($\int_V q \, dV = \eta Q_{\text{tot}}$).

2. **Multi-Phase Metallurgy with `ABQ_PHASE_TRANS`:**
Coupling non-linear thermal material properties (16MnCr5) with Abaqus's built-in phase transformation tables:
- Latent heat across the melting range ($T_{\text{solidus}} \to T_{\text{liquidus}}$).
- Diffusional solid-state transformations (Austenite $\to$ Ferrite/Bainite) via Johnson-Mehl-Avrami (JMA) kinetics.
- Martensitic transformation via Koistinen-Marburger.
- RLS (Raw-Liquid-Solid) state tracking initialized for wrought/machined components.

3. **Time-Stepping Strategy for BTR (Brittleness Temperature Range):**
To accurately capture hot-cracking risk during the high-susceptibility solidification interval, the welding step enforces a fine, uniform time increment ($\Delta t = 0.012\text{ s}$) that extends beyond beam shut-off, before transitioning to adaptive time-stepping (`DELTMX = 25 °C`) during bulk cooling.

4. **Solver-Free Python Verification Suite:**
A Python pre-processing pipeline generates the graded hexahedral mesh (DC3D8) and verifies element Jacobians and Gauss quadrature energy integration ($2\times 2\times 2$ rule) directly on the finite element mesh, proving $<1.9\%$ discretization error prior to solver submission.

Check out the full documentation, equations, and code on GitHub: [link to your repo]

#FEA #Abaqus #Simulation #LaserWelding #Metallurgy #ComputationalMechanics #Fortran #Python
