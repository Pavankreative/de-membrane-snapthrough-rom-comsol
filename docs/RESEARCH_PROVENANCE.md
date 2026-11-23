# Research and software provenance

Pavan Kumar Karsinapu is the manuscript author. The manuscript reports that he built and ran the COMSOL models. The files under `spectral_limits/` include the exported axisymmetric eigenvalues and their numerical post-processing. The shell crossings are stable-side extrapolations. The non-axisymmetric crossing values in Table 7 are predictions from the fitted plate model.

Claude (Anthropic) assisted with manuscript drafting and editing, analysis and post-processing scripts, and COMSOL setup guidance. For revision 17, ChatGPT (OpenAI) implemented and executed the manifold, basin, Lyapunov and reference-trajectory calculations, regenerated figures, compared the computed values with the manuscript, and assisted with interpretation and revision. The revision 17 workflow used the archived COMSOL results. The new dynamics integrations were performed with Python and C.

The [dynamics supplement](../planar_dynamics/README.md) lists the executed programs, inputs, reference environment and archived outputs. Its checksum manifest supports a byte-level comparison with the reference run. The [spectral supplement](../spectral_limits/README.md) identifies the finite-element exports used in each fit and the limits of their interpretation.
