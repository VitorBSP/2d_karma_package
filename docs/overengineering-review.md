# Revisão de utilidade e *overengineering*

Avaliação, componente por componente, do que o trabalho entrega, do que é
útil e do que é excesso. Ela inclui o que eu mesmo construí em excesso durante
a auditoria. Tamanhos medidos no commit atual da branch
`claude/tender-pascal-lmicq4`.

## Veredito

**O pacote em si não está superdimensionado:** são 1.827 linhas de R, uma
única dependência (`stats`) e uma API que cobre o que o artigo usa. O excesso
está **em volta dele**: variantes de pipeline, a maquinaria de validação,
documentação sobreposta, cinco vignettes e 8 mil linhas de código legado.
Esse material se acumulou porque a auditoria precisou comparar várias versões
do mesmo procedimento. Com as perguntas respondidas, cerca de **metade do que
está fora de `R/` pode sair ou ser fundida** sem perder nada que alguém vá
usar.

| Componente | Tamanho | Útil? | Recomendação |
|---|---|---|---|
| Pacote `R/` | 1.827 linhas, 22 exports | Sim | Manter; reduzir exports para 16 (Sec. 1) |
| Testes `tests/` | 980 linhas, 140 testes | Sim | Manter |
| Pipelines `analysis/` (sem legado) | 1.288 linhas | Sim, com variantes demais | Remover 3 variantes (Sec. 2) |
| `validate_results.R` | 270 linhas | Parcialmente | Simplificar (Sec. 3) |
| Código legado `analysis/legacy/` | 8.023 linhas | Só como registro | Mover para uma tag Git (Sec. 5) |
| Vignettes | 5, 1.019 linhas | Sim, mas sobrepostas | Reduzir para 3 (Sec. 4) |
| `docs/` | 4 documentos, 806 linhas | Sim para os autores | Fundir em 2 (Sec. 4) |
| CI | 6 workflows | 4 úteis | Remover 1, rodar 1 (Sec. 3) |
| Dados brutos | 33 rasters NDVI, 1 CSV GRDI | 1 raster usado | Remover 32 rasters (Sec. 5) |

## 1. Pacote (`R/`): o que exportar

O código interno é enxuto e cada parte tem razão de existir: o tratamento das
singularidades removíveis em δ = 1 e δ = 2, `log1p`/`expm1` na distribuição,
a semente local, o escalonamento do otimizador. O excesso está no número de
funções públicas. Cada export é documentação que o CRAN audita e interface
que vocês terão de manter compatível.

| Função | Uso real | Veredito |
|---|---|---|
| `dkuma`, `pkuma`, `qkuma`, `rkuma` | Base de tudo | **Manter** |
| `karma2d_fit`, `karma2d_control`, métodos S3 | Uso principal | **Manter** |
| `karma2d_simulate` | Monte Carlo, exemplos | **Manter** |
| `karma2d_residuals`, `karma2d_detect_anomalies`, `karma2d_metrics` | Aplicações | **Manter** |
| `karma2d_lags` | Documenta a ordem dos coeficientes | **Manter** |
| `karma2d_start` | Necessário para várias partidas em imagens pequenas (máximos locais) | **Manter** |
| `karma2d_loglik`, `karma2d_score`, `karma2d_information` | Verificação por revisores; base dos testes de derivadas | **Manter**, marcadas como "avançadas" |
| `karma2d_diagnostics` | Pequena, útil | **Manter** |
| `karma2d_structuring_element` | Necessária para escolher 3×3 ou 3×5 | **Manter** |
| `karma2d_fisher` | Repete `karma2d_information()` aplicada aos dados do ajuste | **Remover** |
| `morph_erode`, `morph_dilate`, `morph_open`, `morph_close` | Só servem à detecção; os pipelines chamam `morph_close(morph_open(...))` à mão porque falta essa opção | **Tornar internas** e acrescentar `morphology = "opening_closing"` a `karma2d_detect_anomalies()` |

Resultado: de 22 para 16 exports, e a detecção do artigo vira uma chamada só.

**Lacunas de usabilidade** encontradas ao escrever o relatório de uso
(`docs/usage-report.md`), em ordem de importância:

* A abertura seguida de fechamento, que é o pós-processamento do artigo, não
  está disponível em `karma2d_detect_anomalies()` (proposta acima).
* O recorte das primeiras `2w` linhas e colunas depois da predição sem borda
  e o recorte da região efetiva para as métricas de calibração são feitos à
  mão (`m[-(1:2), -(1:2)]`), o que é fácil de errar. Um argumento
  `crop = TRUE` em `predict()` resolveria.
* `predict()` e `karma2d_residuals()` param com qualquer `NA`; em rasters
  reais (GRDI) é preciso imputar antes. Basta documentar a imputação pela
  média usada no artigo, sem criar mais uma opção.
* `padding = "duplicate"` muda também as predições internas quando há parte
  MA; isso agora está documentado em `?predict.karma2d_fit`.

**Opções a decidir com os autores:**

* `predict(padding = "duplicate")` implementa a borda descrita na Sec. 3.1 do
  artigo. A validação do código dos autores mostrou que **essa borda nunca foi
  usada** para gerar os resultados. Se o texto do artigo for corrigido, a
  opção sai: cerca de 25 linhas e 8 testes.
* `karma2d_control(shape_transform = "identity")` só existe para reproduzir o
  caminho do otimizador original. Ela tem um uso legítimo: a segunda partida
  em imagens pequenas. Manter.
* As ligações probit, cloglog e cauchit custam pouco, são testadas e o artigo
  as cita. Manter.

## 2. Pipelines (`analysis/`)

Úteis e enxutos: `run_monte_carlo.R` (fluxos independentes, *checkpoints*),
`summarize_monte_carlo.R` e o benchmark gaussiano.

Excesso acumulado durante a auditoria:

| Item | Por que existe | Veredito |
|---|---|---|
| `procedure = "padded"` em `application.R` | Era o meu primeiro pipeline, antes de chegarem os scripts dos autores | **Remover**; o procedimento do artigo é o único que interessa |
| Linhas gaussianas "como publicado" | Provam a origem dos números da Tabela 5 | **Manter até o artigo ser corrigido**; depois remover |
| Elemento losango 5×5 do EBImage | Hipótese sobre as Figuras 7 e 10, já descartada (os autores usam `mmand` com caixa 3×3) | **Remover** (sobrou no ramo "padded") |
| `--scenario_set=legacy_scripts` no Monte Carlo | Parâmetros dos scripts que **não** geraram as tabelas | **Remover** |
| `legacy_information()` | Reproduz os erros-padrão publicados | **Manter enquanto o artigo usar a matriz mista** |
| Seleção de modelos com ordens 3 | Segue o protocolo do artigo, mas são ajustes de 23 a 32 parâmetros, instáveis entre plataformas, que nunca são escolhidos | Restringir a ordens ≤ 2 (economiza ~30 s por execução) |

Ganho estimado: cerca de 150 linhas a menos em `application.R` e
`config.R`, e um ramo de lógica a menos em toda a validação.

## 3. Validação e CI

**`validate_results.R`** foi útil: detectou corretamente que a mudança do
otimizador alterava os resultados. Mas o custo foi alto: **as quatro falhas do
CI que ele produziu foram alarmes falsos** de diferença numérica entre versões
do R (2e-5 numa seleção de modelos, 1,2e-7 num viés, 3e-5 numa assimetria). Cada uma exigiu um
remendo: tolerância por tabela, filtro de linhas, escala do viés. Os achados
científicos reais (os defeitos do gaussiano, a matriz mista) vieram de análise
manual, não dele.

Simplificação recomendada:

* Manter a checagem de regressão do NDVI (tabelas 4 e 5, detecção), que é
  determinística e barata.
* Comparar o Monte Carlo rápido (20 réplicas) só contra o snapshot, sem a
  checagem "contra o artigo", que nunca reprova.
* Uma tolerância única, documentada, em vez de três mecanismos.
* Meta: cerca de 120 linhas em vez de 270.

**Workflows:**

| Workflow | Veredito |
|---|---|
| R CMD check (5 plataformas) | **Essencial** |
| lint | Útil; falhou uma vez por mudança de versão do lintr |
| validate-results | Útil depois de simplificado |
| pkgdown | Útil quando a branch for para a `main` |
| test-coverage | **Baixo valor**: sem token do Codecov só imprime um número. Remover ou configurar |
| reproduce-paper (manual) | **Nunca foi executado.** Rodar uma vez em modo `quick` ou remover; workflow não testado é promessa falsa |

## 4. Documentação

* **Vignettes, de 5 para 3.** Juntar "simulation-and-inference" em
  "getting-started"; manter "model-theory" e "anomaly-detection"; tirar
  "reproducing-the-paper" do pacote e deixá-la em `analysis/README.md`, porque
  serve aos autores e não ao usuário. Cerca de 360 linhas a menos, e menos
  tempo de compilação no CRAN.
* **`docs/`, de 4 para 2 documentos.** `technical-audit.md` virou um
  palimpsesto, com notas de "atualização" e seções reescritas depois que os
  scripts dos autores chegaram. Recomendo fundi-lo com
  `validation-authors-code.md` num único "relatório para os autores", com as
  perguntas em aberto no topo. Manter `data-provenance.md`;
  `monte-carlo-r1000.md` vira um anexo do relatório.

## 5. Dados e legado

* **32 dos 33 rasters NDVI nunca são usados** (só `NDVI_2019_10.tif`).
  Remover libera 1,5 MB do checkout (o histórico do Git continua com eles).
* **O CSV do GRDI (6,7 MB)** é uma extração do arquivo usado pelo artigo,
  mas **não contém a janela do artigo**: uma busca exaustiva nas 362.604
  posições possíveis encontrou, no melhor caso, um ajuste a 4 erros-padrão
  das estimativas publicadas. Hoje ele serve só para o teste funcional do
  pipeline GRDI. Manter, mas não acrescentar outras cópias.
* **`analysis/legacy/` (8.023 linhas)** tem valor de registro, não de uso. O
  baseline dos testes já está congelado em `tests/testthat/fixtures/`.
  Recomendo criar uma tag `legacy-original` e remover a pasta da `main`;
  `data-raw/legacy-baseline.R` passaria a apontar para a tag.

## 6. O que **não** simplificar

* **Os testes de derivadas e da informação contra integração numérica.** Foram
  eles que encontraram os erros de sinal do Apêndice B.
* **O baseline contra o código original.** É a prova de que a refatoração não
  mudou o modelo.
* **As opções `information = "expected"/"observed"`.** A escolha da matriz é
  uma pergunta em aberto para o artigo.
* **Os fluxos aleatórios independentes e os *checkpoints* do Monte Carlo.**
  Sem eles a execução de 10.000 réplicas não é reprodutível nem retomável.

## 7. Plano proposto

Numa PR separada, depois das decisões dos autores:

1. API: remover `karma2d_fisher`; tornar internas as funções `morph_*`;
   adicionar `morphology = "opening_closing"`.
2. Pipelines: remover `procedure = "padded"`, o losango e
   `legacy_scripts`; limitar a seleção a ordens ≤ 2.
3. Validação: reescrever `validate_results.R` com uma tolerância e sem a
   checagem de Monte Carlo contra o artigo.
4. CI: remover test-coverage (ou configurar Codecov); executar
   reproduce-paper uma vez.
5. Documentação: 3 vignettes; fundir os relatórios.
6. Dados e legado: remover os 32 rasters; mover o legado para uma tag.

Efeito estimado: cerca de 9 mil linhas e 1,5 MB a menos na `main`, 6 exports
a menos, nenhuma perda de funcionalidade usada pelo artigo.
