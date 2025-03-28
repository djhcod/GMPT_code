# User interface ------------------------------------

ui <- page_navbar(
  theme = bs_theme(bootswatch = "morph"),
  inverse = FALSE,
  fillable_mobile = TRUE,
  fluid = FALSE,
  title = "GMPT",


  # 介绍页 ---------------------------------------------------------------------

  nav_panel(
    div(
      img(
        src = "banner_intro.jpg",
        width = "100%",
        style = "display: block"
      ),
      style = "text-align: center;"
    ),
    title = "Intro",
    fluidRow(
      layout_column_wrap(
        card(
          card_header("About this website", style = "font-size: 40px; font-weight: bold"),
          p(
            strong(
              "Gynecologic Malignancies Prediction Toolkit (GMPT)",
              style = "color: rgb(230, 154, 156)"
            ),
            "is a one-stop web application that provides malignancy risk assessment,
        prognosis prediction, and treatment decision-making for gynecologic malignancies.
        The current GMPT encompasses the uterine sarcoma prognostic prediction model,
          the uterine sarcoma adjuvant radiotherapy decision tool,
          as well as the inflammatory marker-based ovarian tumor malignancy risk
          assessment tool and the ovarian cancer prognostic prediction model
          based on the CCL3L1+ Inflam-TAMs-related prognostic gene signature (MRPGS)."
          ),
          p(
            "The development of GMPT was based on large-scale population data,
          clinical case analysis, and single-cell sequencing. We investigated
          the impact of peritoneal cytology, adjuvant radiotherapy, and inflammatory
          response markers on the prognosis of uterine sarcomas and ovarian cancers
          in a multidimensional manner and constructed a series of highly efficient
          predictive models. We revealed the independent prognostic value of
          malignant peritoneal cytology for uterine sarcoma and developed a
          risk stratification system that can guide radiotherapy decisions.
          The inflammatory marker-based column-line diagram model and the
          CCL3L1+ Inflam-TAMs-related gene-based integrated prognostic model
          provided new options for malignant risk assessment and prognostic
          prediction of ovarian cancer. Finally, the GMPT web application that
          integrates all models achieves clinical translation, and its
          open-source methodology provides a standardized tool for
          individualized treatment strategy development and prognostic
          assessment."
          ),
          p(
            "The GMPT web application will be maintained long-term. In the future,
          we will continue enhancing the GMPT website to include a broader range
          of gynecologic malignancies and predictive tools, aiming to provide a
          comprehensive, one-stop prediction platform."
          )
        )
      )
    ),
    fluidRow(
      column(3, img(src = "tree.png", style = "width: 100%; padding: 10px;")),
      column(3, img(src = "lasso.png", style = "width: 100%; padding: 10px;")),
      column(3, img(src = "rsf.png", style = "width: 100%; padding: 10px;")),
      column(3, img(src = "sc.png", style = "width: 100%; padding: 10px;"))
    ),
    fluidRow(
      layout_column_wrap(
        card(
          card_header("Links", style = "font-size: 30px; font-weight: bold"),
          p(
            icon("github", class = "fa-brands fa-github", lib = "font-awesome", style = "color: rgb(230, 154, 156)"),
            strong("Source code:  ", style = "color: rgb(230, 154, 156)"),
            a(
              "https://github.com/djhcod/GMPT_code",
              href = "https://github.com/djhcod/GMPT_code",
              target = "_blank",
            ),
          ),
          p(
            icon("app-store", class = "fa-brands fa-app-store", lib = "font-awesome", style = "color: rgb(230, 154, 156)"),
            strong("GMPT desktop application:  ", style = "color: rgb(230, 154, 156)"),
            a(
              "https://github.com/djhcod/GMPT_code/releases",
              href = "https://github.com/djhcod/GMPT_code/releases",
              target = "_blank",
            ),
          ),
          p(
            icon("laptop-code", class = "fa-solid fa-laptop-code", lib = "font-awesome", style = "color: rgb(230, 154, 156)"),
            strong("Local deployment of GMPT:  ", style = "color: rgb(230, 154, 156)"),
            a(
              "https://github.com/djhcod/GMPT_code/tree/main/web_nomogram",
              href = "https://github.com/djhcod/GMPT_code/tree/main/web_nomogram",
              target = "_blank",
            ),
            "(run through RStudio after downloading all files completely)"
          )
        )
      )
    ),
    fluidRow(
      layout_column_wrap(
        card(
          card_header("SCI publications to date for this project:", style = "font-size: 20px; font-weight: bold"),
          p(
            "Du J, Hu D, Xing Y, et al. Correlation between malignant peritoneal cytology
            and survival in patients with uterine leiomyosarcoma and endometrial
            stromal sarcoma[J]. International Journal of Gynecologic Cancer,
            2024, 34(1): 58-65. DOI: ",
            a(
              "10.1136/ijgc-2023-004792.",
              href = "https://doi.org/10.1136/ijgc-2023-004792",
              target = "_blank",
            )
          ),
          p(
            "Du J, Cheng Y, Hu D, et al. A nomogram-based overall survival
            stratification to identify uterine sarcoma patients without
            distant metastases who may benefit from adjuvant radiotherapy[J].
            Gynecologic Oncology, 2023, 169: 17-26. DOI: ",
            a(
              "10.1016/j.ygyno.2022.11.023.",
              href = "https://doi.org/10.1016/j.ygyno.2022.11.023",
              target = "_blank",
            )
          )
        )
      )
    ),
    p("Du Junhong Copyright 2025",
      style = "text-align: center"
    )
  ),



  # page 1 ------------------------------------------------------------------

  nav_panel(
    title = "Peritoneal cytology-based model",
    h1(
      style = "color: black;
              text-align: left;
              background-image: url('banner_peri.jpg');
              background-size: cover;
              padding: 100px",
      strong("Peritoneal cytology"),
      "-based overall survival prediction model for patients with uterine sarcoma",
      br(),
      a(
        href = "https://github.com/djhcod/GMPT_code/tree/main/sarcoma_peri",
        icon("github"),
        style = "color: black",
        target = "_blank"
      ),
      a(
        href = "http://dx.doi.org/10.1136/ijgc-2023-004792",
        icon("paperclip"),
        style = "color: black",
        target = "_blank"
      )
    ),
    fluidRow(
      layout_column_wrap(
        card(
          card_header("Introduction", style = "font-size: 28px; font-weight: bold"),
          p(
            "Our ",
            a(
              "previous study",
              href = "http://dx.doi.org/10.1136/ijgc-2023-004792",
              target = "_blank",
            ),
            " confirmed that malignant peritoneal cytology correlated with a
      reduced 5-year overall survival (OS) rate and median survival time in patients with
      uterine leiomyosarcoma and endometrial stromal sarcoma. Therefore, on the basis of
      previous studies we constructed this web calculator based on peritoneal
      cytologic status as well as other significant clinical parameters for
      predicting overall survival in patients with uterine leiomyosarcoma and
      endometrial stromal sarcoma."
          ),
          p(
            icon("circle-info", class = "fa-solid fa-circle-info", lib = "font-awesome"),
            "NOTE: only for patients with uterine leiomyosarcoma and endometrial stromal sarcoma",
            style = "color: rgb(230, 154, 156); font-weight: bold"
          )
        )
      )
    ),


    # 输入 -------------------------------------------------------

    h1("Input", style = "font-size: 28px; font-weight: bold"),
    fluidRow(
      layout_column_wrap(
        card(
          card_header("Age (years)", style = "font-size: 20px; font-weight: bold"),
          sliderInput(
            inputId = "age_peri",
            label = NULL,
            min = 18,
            max = 90,
            value = 50
          )
        ),
        card(
          card_header("Grade", style = "font-size: 20px; font-weight: bold"),
          radioButtons(
            inputId = "grade2_peri",
            label = NULL,
            choices = list(
              "Low-grade" = 1,
              "High-grade" = 2,
              "Gx" = 3
            )
          )
        ),
        card(
          card_header("FIGO stage", style = "font-size: 20px; font-weight: bold"),
          radioButtons(
            inputId = "figo_peri",
            label = NULL,
            choices = c(
              "Stage I" = 1,
              "Stage II" = 2,
              "Stage III" = 3,
              "Stage IV" = 4
            )
          )
        )
      )
    ),
    fluidRow(
      layout_column_wrap(
        card(
          card_header("Adjuvant chemotherapy", style = "font-size: 20px; font-weight: bold"),
          radioButtons(
            inputId = "chem_peri",
            label = NULL,
            choices = c(
              "No" = 0,
              "Yes" = 1
            )
          )
        ),
        card(
          card_header("Peritoneal cytology status", style = "font-size: 20px; font-weight: bold"),
          radioButtons(
            inputId = "peri",
            label = NULL,
            choices = c(
              "Negtive" = 0,
              "Malignant" = 1
            )
          )
        ),
        card(
          div(
            style = "
            display: flex;
            justify-content: center;
            align-items: center;
            height: 100%
            ",
            actionButton(
              inputId = "action_peri",
              label = "Predict!",
              icon("paper-plane", class = "fa-solid fa-paper-plane", lib = "font-awesome"),
              style = "
              font-weight: bold;
              background-color: rgb(230, 154, 156);
              border-color: rgb(200, 120, 120);
              color: white;
              box-shadow: 0px 4px 10px rgba(230, 154, 156, 0.7)
              "
            )
          )
        )
      )
    ),


    # 输出 --------------------------------------------------

    h1("Output", style = "font-size: 28px; font-weight: bold"),

    # Case report 独占一行
    fluidRow(
      column(
        width = 12,
        card(
          card_title("Case report"),
          tableOutput(outputId = "case_report_peri")
        )
      )
    ),

    # 下一行，左侧 "Predicted OS curve"，右侧 "Total points" 和 "Predicted OS"
    fluidRow(
      column(
        width = 7,
        card(
          card_title("Predicted OS curve"),
          plotOutput(outputId = "os_curve_peri"),
          style = "height: 100%"
        )
      ),
      column(
        width = 5,
        div(
          style = "display: flex; flex-direction: column; height: 100%",
          card(
            card_title("Total points"),
            div(
              style = "color: rgb(230, 154, 156); font-size: 40px; font-weight: bold",
              tableOutput(outputId = "points_peri")
            ),
            style = "flex: 1"
          ),
          card(
            card_title("Predicted OS"),
            tableOutput(outputId = "predictions_peri"),
            style = "flex: 2"
          )
        )
      )
    ),
    p("Du Junhong Copyright. 2025", style = "text-align: center")
  ),





  # page 2 ------------------------------------------------------------------

  nav_panel(
    title = "Radiotherapy Decision",
    h1(
      "A nomogram-based overall survival stratification to identify uterine sarcoma
       patients without distant metastases",
      strong("who may benefit from adjuvant radiotherapy"),
      br(),
      a(
        href = "https://github.com/djhcod/GMPT_code/tree/main/sarcoma_rad",
        icon("github"),
        style = "color: black",
        target = "_blank"
        ),
      a(
        href = "https://doi.org/10.1016/j.ygyno.2022.11.023",
        icon("paperclip"),
        style = "color: black",
        target = "_blank"
        ),
      style = "color: black;
              text-align: left;
              background-image: url('banner_rad.jpeg');
              background-size: cover;
              padding: 100px"
    ),
    fluidRow(
      layout_column_wrap(
        card(
          card_header("Introduction", style = "font-size: 28px; font-weight: bold"),
          p(
            "Adjuvant radiotherapy has been commonly performed in uterine sarcoma
            patients, but its role in overall survival (OS) remains controversial.
            Therefore, we build this online prognostic stratification tool to
            identify uterine sarcoma patients who might benefit from adjuvant radiotherapy.
            The tool is able to predict the overall survival (OS) of a given
            uterine sarcoma patient in the absence of radiotherapy using a
            few simple clinical variables, and categorizes the patient into low,
            intermediate, or high risk groups based on model scores."
          ),
          p(
            "According to our ",
            a(
              "published research",
              href = "https://doi.org/10.1016/j.ygyno.2022.11.023",
              target = "_blank",
            ),
            "for patients in the high-risk group, adjuvant radiotherapy significantly
            improved the 5-year OS and median survival time by 26.4% and 17 months,
            respectively (P < 0.001); while radiotherapy failed to improve the
            survival outcomes of patients in the low- and intermediate-risk groups."
          ),
          p(
            icon("circle-info", class = "fa-solid fa-circle-info", lib = "font-awesome"),
            "NOTE: only for uterine sarcoma patients without distant metastases",
            style = "color: rgb(230, 154, 156); font-weight: bold"
          )
        )
      )
    ),
    h1("Input", style = "font-size: 28px; font-weight: bold"),
    fluidRow(
      layout_column_wrap(
        card(
          card_header("Age (years)", style = "font-size: 20px; font-weight: bold"),
          sliderInput(
            inputId = "age",
            label = NULL,
            min = 18,
            max = 90,
            value = 50
          )
        ),
        card(
          card_header("Grade", style = "font-size: 20px; font-weight: bold"),
          radioButtons(
            inputId = "grade",
            label = NULL,
            choices = list(
              "Low-grade" = 1,
              "High-grade" = 2,
              "Gx" = 3
            )
          )
        ),
        card(
          card_header("Tumor size", style = "font-size: 20px; font-weight: bold"),
          radioButtons(
            inputId = "tumor_size",
            label = NULL,
            choices = list(
              "<70 mm" = 1,
              "70–165 mm" = 2,
              ">165 mm" = 3,
              "Unspecific" = 4
            )
          )
        ),
        card(
          card_header("Histology", style = "font-size: 20px; font-weight: bold"),
          radioButtons(
            inputId = "his",
            label = NULL,
            choices = list(
              "Leiomyosarcoma" = 1,
              "Endometrial stromal sarcoma" = 2,
              "Adenosarcoma" = 3
            )
          )
        )
      )
    ),
    fluidRow(
      layout_column_wrap(
        card(
          card_header("T stage", style = "font-size: 20px; font-weight: bold"),
          radioButtons(
            inputId = "T_stage",
            label = NULL,
            choices = c(
              "T1" = 1,
              "T2" = 2,
              "T3" = 3,
              "T4" = 4
            )
          )
        ),
        card(
          card_header("N stage", style = "font-size: 20px; font-weight: bold"),
          radioButtons(
            inputId = "N_stage",
            label = NULL,
            choices = c(
              "N0" = 0,
              "N1" = 1
            )
          )
        ),
        card(
          card_header("Adjuvant chemotherapy", style = "font-size: 20px; font-weight: bold"),
          radioButtons(
            inputId = "chemotherapy",
            label = NULL,
            choices = c(
              "No" = 0,
              "Yes" = 1
            )
          )
        ),
        card(
          div(
            style = "
            display: flex;
            justify-content: center;
            align-items: center;
            height: 100%
            ",
            actionButton(
              inputId = "action_rad",
              label = "Predict!",
              icon("paper-plane", class = "fa-solid fa-paper-plane", lib = "font-awesome"),
              style = "
              font-weight: bold;
              background-color: rgb(230, 154, 156);
              border-color: rgb(200, 120, 120);
              color: white;
              box-shadow: 0px 4px 10px rgba(230, 154, 156, 0.7)
              "
            )
          )
        )
      )
    ),
    h1("Output", style = "font-size: 28px; font-weight: bold"),
    fluidRow(
      column(
        width = 12,
        card(
          card_title("Case report"),
          tableOutput(outputId = "case_report_rad")
        )
      )
    ),
    fluidRow(
      column(
        width = 7,
        card(
          card_title("Predicted OS curve"),
          plotOutput(outputId = "os_curve_rad"),
          style = "flex: 10"
        ),
        card(
          card_title("Predicted OS without radiotherapy"),
          tableOutput(outputId = "predictions_rad"),
          style = "flex: 10"
        )
      ),
      column(
        width = 5,
        div(
          style = "display: flex; flex-direction: column; height: 100%",
          card(
            card_title("Total points"),
            div(
              style = "color: rgb(230, 154, 156); font-size: 40px; font-weight: bold",
              textOutput(outputId = "points_rad")
            ),
            style = "flex: 5"
          ),
          card(
            card_title("Risk group"),
            div(
              style = "color: rgb(230, 154, 156); font-size: 30px; font-weight: bold",
              uiOutput(outputId = "risk_rad")
            ),
            style = "flex: 5"
          ),
          card(
            card_title("Adjuvant benefit"),
            div(
              style = "color: rgb(230, 154, 156); font-size: 30px; font-weight: bold",
              uiOutput(outputId = "benefit_rad")
            ),
            style = "flex: 10"
          )
        )
      )
    ),
    p(
      icon("circle-info", class = "fa-solid fa-circle-info", lib = "font-awesome"),
      "This tool aids in decision-making for adjuvant radiotherapy in
            uterine sarcoma patients; however, the final decision should rely on
            the oncologist's clinical assessment of the patient's general conditions,
            comorbidities, intraoperative condition, and socioeconomic factors.",
      style = "text-align: center"
    ),
    hr(),
    p("Du Junhong Copyright 2025",
      style = "text-align: center"
    )
  ),









  # page 3 ------------------------------------------------------------------

  nav_panel(
    title = "OC Predict",
    h1(
      style = "color: black;
              text-align: left;
              background-image: url('banner_oc.jpg');
              background-size: cover;
              padding: 100px",
      strong("Systemic Inflammatory Response Markers"),
      "-based Nomogram to Estimate the Malignant Probability in Patients with Ovarian Masses",
      br(),
      a(
        href = "https://github.com/djhcod/GMPT_code/tree/main/inflammatory_markers",
        icon("github"),
        style = "color: black",
        target = "_blank"
      )
    ),
    fluidRow(
      layout_column_wrap(
        card(
          card_header("Introduction", style = "font-size: 28px; font-weight: bold"),
          p(
            "The model study included 1087 patients who underwent surgical treatment
            for ovarian tumors. Firstly, systemic inflammatory response markers
            that can effectively differentiate between benign and malignant ovarian
            masses were identified based on LASSO-logistic regression.
            Three levels of models were established based on different variables (Models 1-3).
            The inclusion of two inflammatory indicators in Model 3, namely fibrinogen (FIB)
            and monocyte-to-lymphocyte ratio (MLR), was demonstrated to significantly improve
            the accuracy of the model in differentiating between benign and malignant ovarian
            masses by calculating ΔAUC, NRI and IDI. Finally, a column-line graph prediction
            model was constructed based on Model 3 and demonstrated to have higher accuracy
            in predicting ovarian malignancy and identifying early ovarian cancer compared
            with RMI, LR2, ROMA, CPH-I and R-OPS models."
          )
        )
      )
    ),


    # 输入 -------------------------------------------------------

    h1("Input", style = "font-size: 28px; font-weight: bold"),
    fluidRow(
      layout_column_wrap(
        card(
          card_header("Age (years)", style = "font-size: 20px; font-weight: bold"),
          sliderInput(
            inputId = "age_ov",
            label = NULL,
            min = 18,
            max = 90,
            value = 50
          )
        ),
        card(
          card_header("CA125 (U/ml)", style = "font-size: 20px; font-weight: bold"),
          radioButtons(
            inputId = "CA125_ov",
            label = NULL,
            choices = c(
              "≤ 100" = 0,
              "> 100" = 1
            )
          )
        ),
        card(
          card_header("HE4 (U/ml)", style = "font-size: 20px; font-weight: bold"),
          radioButtons(
            inputId = "HE4_ov",
            label = NULL,
            choices = c(
              "≤ 64.5" = 0,
              "> 64.5" = 1
            )
          )
        ),
        card(
          card_header("Fibrinogen (g/L)", style = "font-size: 20px; font-weight: bold"),
          radioButtons(
            inputId = "FIB_ov",
            label = NULL,
            choices = c(
              "≤ 3.28" = 0,
              "> 3.28" = 1
            )
          )
        )
      )
    ),
    fluidRow(
      layout_column_wrap(
        card(
          card_header("Monocyte (10^9/L)", style = "font-size: 20px; font-weight: bold"),
          sliderInput(
            inputId = "mon_ov",
            label = NULL,
            min = 0.1,
            max = 2,
            value = 0.32
          )
        ),
        card(
          card_header("Lymphocyte (10^9/L)", style = "font-size: 20px; font-weight: bold"),
          sliderInput(
            inputId = "lym_ov",
            label = NULL,
            min = 0.1,
            max = 3,
            value = 1.48
          )
        ),
        card(
          card_header("Strong blood flow signal", style = "font-size: 20px; font-weight: bold"),
          radioButtons(
            inputId = "Bloodflow_ov",
            label = NULL,
            choices = c(
              "No" = 0,
              "Yes" = 1
            )
          )
        ),
        card(
          card_header("Presence of solid areas", style = "font-size: 20px; font-weight: bold"),
          radioButtons(
            inputId = "Solid_areas_ov",
            label = NULL,
            choices = c(
              "No" = 0,
              "Yes" = 1
            )
          )
        )
      )
    ),
    fluidRow(
      layout_column_wrap(
        div(
          style = "
            display: flex;
            justify-content: center;
            align-items: center;
            height: 100%
            ",
          actionButton(
            inputId = "action_ov",
            label = "Predict!",
            icon("paper-plane", class = "fa-solid fa-paper-plane", lib = "font-awesome"),
            style = "
              font-weight: bold;
              background-color: rgb(230, 154, 156);
              border-color: rgb(200, 120, 120);
              color: white;
              box-shadow: 0px 4px 10px rgba(230, 154, 156, 0.7)
              "
          )
        )
      )
    ),







    # 输出 --------------------------------------------------

    h1("Output", style = "font-size: 28px; font-weight: bold"),

    # Case report 独占一行
    fluidRow(
      column(
        width = 12,
        card(
          card_title("Case report"),
          tableOutput(outputId = "case_report_ov")
        )
      )
    ),
    fluidRow(
      column(
        width = 6,
        card(
          card_title("Total points"),
          div(
            style = "color: rgb(230, 154, 156); font-size: 40px; font-weight: bold",
            textOutput(outputId = "points_ov")
          )
        )
      ),
      column(
        width = 6,
        card(
          card_title("Malignant probability"),
          div(
            style = "color: rgb(230, 154, 156); font-size: 40px; font-weight: bold",
            uiOutput(outputId = "predict_ov")
          )
        )
      )
    ),
    p("Du Junhong Copyright. 2025", style = "text-align: center")
  ),








  # page 4 ------------------------------------------------------------------

  nav_panel(
    title = "MRPGS",
    h1(
      style = "color: black;
              text-align: left;
              background-image: url('banner_sc.jpg');
              background-size: cover;
              padding: 100px",
      strong("CCL3L1+ Inflam-TAMs"),
      " related prognostic gene signature (MRPGS) -based prognostic prediction
      model for ovarian cancer.",
      br(),
      a(
        href = "https://github.com/djhcod/GMPT_code/tree/main/macrophage_oc",
        icon("github"),
        style = "color: black",
        target = "_blank"
      )
    ),
    fluidRow(
      layout_column_wrap(
        card(
          card_header("Introduction", style = "font-size: 28px; font-weight: bold"),
          p(
            "Our ",
            a(
              "previous study",
              href = "http://dx.doi.org/10.1136/ijgc-2023-004792",
              target = "_blank",
            ),
            " identified several novel molecular subtypes of macrophages by
            targeting the monocyte-macrophage system in the MLR with single-cell
            sequencing data and identified CCL3L1+ Inflam-TAMs as a key
            macrophage subpopulation that contributes significantly to
            the prognostic prediction of ovarian cancer, and prognostic
            scoring and comprehensive modeling based on
            CCL3L1+ Inflam-TAMs-related genes showed high accuracy in prognostic
            prediction of ovarian cancer patients."
          ),
        )
      )
    ),


    # 输入 -------------------------------------------------------

    h1("Input", style = "font-size: 28px; font-weight: bold"),
    h2("Gene expression level input", style = "font-size: 20px; font-weight: bold"),
    fluidRow(
      layout_column_wrap(
        card(
          card_header("C5AR1", style = "font-size: 20px; font-weight: bold"),
          numericInput(
            inputId = "C5AR1_sc",
            label = NULL,
            value = 1000,
            min = 0
          )
        ),
        card(
          card_header("EREG", style = "font-size: 20px; font-weight: bold"),
          numericInput(
            inputId = "EREG_sc",
            label = NULL,
            value = 1000,
            min = 0
          )
        ),
        card(
          card_header("FCGR2A", style = "font-size: 20px; font-weight: bold"),
          numericInput(
            inputId = "FCGR2A_sc",
            label = NULL,
            value = 1000,
            min = 0
          )
        ),
        card(
          card_header("PLK3", style = "font-size: 20px; font-weight: bold"),
          numericInput(
            inputId = "PLK3_sc",
            label = NULL,
            value = 1000,
            min = 0
          )
        ),
        card(
          card_header("SAT1", style = "font-size: 20px; font-weight: bold"),
          numericInput(
            inputId = "SAT1_sc",
            label = NULL,
            value = 1000,
            min = 0
          )
        )
      )
    ),
    fluidRow(
      layout_column_wrap(
        card(
          card_header("SLC11A2", style = "font-size: 20px; font-weight: bold"),
          numericInput(
            inputId = "SLC11A2_sc",
            label = NULL,
            value = 1000,
            min = 0
          )
        ),
        card(
          card_header("SPAG9", style = "font-size: 20px; font-weight: bold"),
          numericInput(
            inputId = "SPAG9_sc",
            label = NULL,
            value = 1000,
            min = 0
          )
        ),
        card(
          card_header("SRGAP1", style = "font-size: 20px; font-weight: bold"),
          numericInput(
            inputId = "SRGAP1_sc",
            label = NULL,
            value = 1000,
            min = 0
          )
        ),
        card(
          card_header("TAOK3", style = "font-size: 20px; font-weight: bold"),
          numericInput(
            inputId = "TAOK3_sc",
            label = NULL,
            value = 1000,
            min = 0
          )
        ),
        card(
          card_header("XBP1", style = "font-size: 20px; font-weight: bold"),
          numericInput(
            inputId = "XBP1_sc",
            label = NULL,
            value = 1000,
            min = 0
          )
        )
      )
    ),
    h2("Clinical variables input", style = "font-size: 20px; font-weight: bold"),
    fluidRow(
      layout_column_wrap(
        card(
          card_header("Age (years)", style = "font-size: 20px; font-weight: bold"),
          sliderInput(
            inputId = "age_sc",
            label = NULL,
            min = 18,
            max = 90,
            value = 50
          )
        ),
        card(
          card_header("Residual disease", style = "font-size: 20px; font-weight: bold"),
          radioButtons(
            inputId = "residual_sc",
            label = NULL,
            choices = list(
              "No Macroscopic disease",
              "1-10 mm",
              "11-20 mm",
              ">20 mm"
            )
          )
        ),
        card(
          card_header("Chemotherapy regimen", style = "font-size: 20px; font-weight: bold"),
          radioButtons(
            inputId = "chemo_sc",
            label = NULL,
            choices = list(
              "Paclitaxel + Carboplatin" = "Paclitaxel + Carboplatin",
              "No chemotherapy" = "No chemotherapy",
              "Paclitaxel + Cisplatin" = "Paclitaxel + Cisplatin",
              "More than 2 drugs" = "More than 2 drugs",
              "Others" = "Others"
            )
          )
        ),
        card(
          card_header("Bevacizumab treatment", style = "font-size: 20px; font-weight: bold"),
          radioButtons(
            inputId = "beva_sc",
            label = NULL,
            choices = list(
              "No" = "No",
              "Yes" = "Yes"
            )
          )
        ),
        card(
          div(
            style = "
            display: flex;
            justify-content: center;
            align-items: center;
            height: 100%
            ",
            actionButton(
              inputId = "action_sc",
              label = "Predict!",
              icon("paper-plane", class = "fa-solid fa-paper-plane", lib = "font-awesome"),
              style = "
              font-weight: bold;
              background-color: rgb(230, 154, 156);
              border-color: rgb(200, 120, 120);
              color: white;
              box-shadow: 0px 4px 10px rgba(230, 154, 156, 0.7)
              "
            )
          )
        )
      )
    ),



    # 输出 --------------------------------------------------

    h1("Output", style = "font-size: 28px; font-weight: bold"),

    # Case report 独占一行
    fluidRow(
      column(
        width = 12,
        card(
          card_title("Case report"),
          tableOutput(outputId = "case_report_sc")
        )
      )
    ),

    # 下一行，左侧 "Predicted OS curve"，右侧 "Total points" 和 "Predicted OS"
    fluidRow(
      column(
        width = 7,
        card(
          card_title("Predicted OS curve"),
          plotOutput(outputId = "os_curve_sc"),
          style = "height: 100%"
        )
      ),
      column(
        width = 5,
        div(
          style = "display: flex; flex-direction: column; height: 100%",
          card(
            card_title("MRPGS"),
            div(
              style = "color: rgb(230, 154, 156); font-size: 40px; font-weight: bold",
              tableOutput(outputId = "points_sc")
            ),
            style = "flex: 1"
          ),
          card(
            card_title("Predicted OS"),
            tableOutput(outputId = "predictions_sc"),
            style = "flex: 2"
          )
        )
      )
    ),
    p("Du Junhong Copyright. 2025", style = "text-align: center")
  ),
  nav_spacer(),
  nav_menu(
    title = "Links",
    align = "right",
    nav_item(tags$a("GitHub repository", href = "https://github.com/djhcod/GMPT_code", target = "_blank")),
    nav_item(tags$a("Desktop application", href = "https://github.com/djhcod/GMPT_code/releases", target = "_blank"))
  )
)
