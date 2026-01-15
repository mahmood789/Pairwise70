"""
Comprehensive Selenium Test for MAFI Calculator Complete Edition
Tests all functions and UI elements in Chrome browser
"""

import time
import json
import os
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import Select
from webdriver_manager.chrome import ChromeDriverManager

# Test results tracking
test_results = {
    "passed": 0,
    "failed": 0,
    "tests": []
}

def log_test(name, passed, details=""):
    """Log test result"""
    status = "PASS" if passed else "FAIL"
    test_results["tests"].append({
        "name": name,
        "status": status,
        "details": details
    })
    if passed:
        test_results["passed"] += 1
    else:
        test_results["failed"] += 1
    print(f"  [{status}] {name}" + (f" - {details}" if details and not passed else ""))

def setup_driver():
    """Setup Chrome driver"""
    options = Options()
    options.add_argument("--start-maximized")
    options.add_argument("--disable-extensions")
    options.add_argument("--disable-popup-blocking")
    options.add_experimental_option("excludeSwitches", ["enable-logging"])

    service = Service(ChromeDriverManager().install())
    driver = webdriver.Chrome(service=service, options=options)
    return driver

def wait_for_element(driver, by, value, timeout=10):
    """Wait for element to be present"""
    return WebDriverWait(driver, timeout).until(
        EC.presence_of_element_located((by, value))
    )

def wait_for_clickable(driver, by, value, timeout=10):
    """Wait for element to be clickable"""
    return WebDriverWait(driver, timeout).until(
        EC.element_to_be_clickable((by, value))
    )

# ============================================================
# TEST FUNCTIONS
# ============================================================

def test_page_load(driver):
    """Test 1: Page loads correctly"""
    print("\n=== TEST GROUP 1: Page Load ===")

    try:
        title = driver.title
        log_test("Page title present", "MAFI Calculator" in title, f"Title: {title}")
    except Exception as e:
        log_test("Page title present", False, str(e))

    try:
        header = driver.find_element(By.TAG_NAME, "header")
        log_test("Header element exists", header is not None)
    except Exception as e:
        log_test("Header element exists", False, str(e))

    try:
        tabs = driver.find_elements(By.CLASS_NAME, "tab")
        log_test("Navigation tabs exist", len(tabs) == 4, f"Found {len(tabs)} tabs")
    except Exception as e:
        log_test("Navigation tabs exist", False, str(e))

def test_navigation_tabs(driver):
    """Test 2: Tab navigation works"""
    print("\n=== TEST GROUP 2: Navigation Tabs ===")

    tab_names = ["calculator", "casestudies", "grade", "documentation"]

    for tab_name in tab_names:
        try:
            # Find and click the tab
            tabs = driver.find_elements(By.CLASS_NAME, "tab")
            target_tab = None
            for tab in tabs:
                onclick = tab.get_attribute("onclick") or ""
                if tab_name in onclick:
                    target_tab = tab
                    break

            if target_tab:
                target_tab.click()
                time.sleep(0.3)

                # Check if corresponding content is visible
                tab_content = driver.find_element(By.ID, f"{tab_name}Tab")
                is_active = "active" in tab_content.get_attribute("class")
                log_test(f"Tab '{tab_name}' navigation", is_active)
            else:
                log_test(f"Tab '{tab_name}' navigation", False, "Tab not found")
        except Exception as e:
            log_test(f"Tab '{tab_name}' navigation", False, str(e))

    # Return to calculator tab
    driver.find_elements(By.CLASS_NAME, "tab")[0].click()
    time.sleep(0.3)

def test_study_management(driver):
    """Test 3: Study management functions"""
    print("\n=== TEST GROUP 3: Study Management ===")

    # Test addStudyRow
    try:
        initial_rows = len(driver.find_elements(By.CSS_SELECTOR, "#studyTableBody tr"))
        add_btn = driver.find_element(By.XPATH, "//button[contains(text(), '+ Add Study')]")
        add_btn.click()
        time.sleep(0.2)
        new_rows = len(driver.find_elements(By.CSS_SELECTOR, "#studyTableBody tr"))
        log_test("addStudyRow() adds one row", new_rows == initial_rows + 1)
    except Exception as e:
        log_test("addStudyRow() adds one row", False, str(e))

    # Test addMultipleRows
    try:
        initial_rows = len(driver.find_elements(By.CSS_SELECTOR, "#studyTableBody tr"))
        add5_btn = driver.find_element(By.XPATH, "//button[contains(text(), '+ Add 5 Studies')]")
        add5_btn.click()
        time.sleep(0.2)
        new_rows = len(driver.find_elements(By.CSS_SELECTOR, "#studyTableBody tr"))
        log_test("addMultipleRows() adds 5 rows", new_rows == initial_rows + 5)
    except Exception as e:
        log_test("addMultipleRows() adds 5 rows", False, str(e))

    # Test clearAllStudies
    try:
        clear_btn = driver.find_element(By.XPATH, "//button[contains(text(), 'Clear All')]")
        clear_btn.click()
        time.sleep(0.2)
        rows = len(driver.find_elements(By.CSS_SELECTOR, "#studyTableBody tr"))
        log_test("clearAllStudies() clears table", rows == 2)  # Should leave 2 empty rows
    except Exception as e:
        log_test("clearAllStudies() clears table", False, str(e))

    # Test remove button
    try:
        remove_btns = driver.find_elements(By.CLASS_NAME, "remove-btn")
        if len(remove_btns) > 0:
            initial_count = len(remove_btns)
            remove_btns[0].click()
            time.sleep(0.2)
            new_count = len(driver.find_elements(By.CLASS_NAME, "remove-btn"))
            log_test("Remove button works", new_count == initial_count - 1)
        else:
            log_test("Remove button works", False, "No remove buttons found")
    except Exception as e:
        log_test("Remove button works", False, str(e))

def test_effect_measure_dropdown(driver):
    """Test 4: Effect measure dropdown"""
    print("\n=== TEST GROUP 4: Effect Measure Dropdown ===")

    try:
        dropdown = Select(driver.find_element(By.ID, "effectMeasure"))
        options = [opt.get_attribute("value") for opt in dropdown.options]
        expected = ["SMD", "MD", "OR", "RR", "HR"]
        log_test("All effect measures present", all(e in options for e in expected), f"Found: {options}")

        # Test selection
        for measure in expected:
            dropdown.select_by_value(measure)
            time.sleep(0.1)
            selected = dropdown.first_selected_option.get_attribute("value")
            log_test(f"Select {measure}", selected == measure)

        # Reset to SMD
        dropdown.select_by_value("SMD")
    except Exception as e:
        log_test("Effect measure dropdown", False, str(e))

def test_clinical_threshold_input(driver):
    """Test 5: Clinical threshold input"""
    print("\n=== TEST GROUP 5: Clinical Threshold Input ===")

    try:
        input_field = driver.find_element(By.ID, "clinicalThreshold")
        input_field.clear()
        input_field.send_keys("0.5")
        value = input_field.get_attribute("value")
        log_test("Clinical threshold accepts value", value == "0.5")

        # Test decimal input
        input_field.clear()
        input_field.send_keys("1.25")
        value = input_field.get_attribute("value")
        log_test("Clinical threshold accepts decimal", value == "1.25")
    except Exception as e:
        log_test("Clinical threshold input", False, str(e))

def test_case_studies(driver):
    """Test 6: Case study loading"""
    print("\n=== TEST GROUP 6: Case Studies ===")

    case_studies = ["antidepressants", "mortality", "pain", "vaccine", "fragile"]

    for cs_name in case_studies:
        try:
            # Go to case studies tab
            tabs = driver.find_elements(By.CLASS_NAME, "tab")
            tabs[1].click()  # Case studies tab
            time.sleep(0.3)

            # Find and click the case study
            case_study_elements = driver.find_elements(By.CLASS_NAME, "case-study")
            found = False
            for cs in case_study_elements:
                onclick = cs.get_attribute("onclick") or ""
                if cs_name in onclick:
                    cs.click()
                    time.sleep(0.3)
                    found = True
                    break

            if not found:
                log_test(f"Load case study: {cs_name}", False, "Case study not found")
                continue

            # Check that studies were loaded
            rows = driver.find_elements(By.CSS_SELECTOR, "#studyTableBody tr")
            has_data = len(rows) >= 3  # At least 3 studies

            # Check first row has data
            if rows:
                first_input = rows[0].find_element(By.CSS_SELECTOR, "input[type='text']")
                has_name = len(first_input.get_attribute("value")) > 0
            else:
                has_name = False

            log_test(f"Load case study: {cs_name}", has_data and has_name, f"{len(rows)} studies loaded")

        except Exception as e:
            log_test(f"Load case study: {cs_name}", False, str(e))

    # Return to calculator tab
    driver.find_elements(By.CLASS_NAME, "tab")[0].click()
    time.sleep(0.3)

def test_load_example_data(driver):
    """Test 7: Load example data button"""
    print("\n=== TEST GROUP 7: Load Example Data ===")

    try:
        # Click Load Example button
        load_btn = driver.find_element(By.XPATH, "//button[contains(text(), 'Load Example')]")
        load_btn.click()
        time.sleep(0.3)

        # Check that data was loaded
        rows = driver.find_elements(By.CSS_SELECTOR, "#studyTableBody tr")
        log_test("Load Example loads data", len(rows) >= 3, f"{len(rows)} studies")

        # Check effect measure was set
        dropdown = Select(driver.find_element(By.ID, "effectMeasure"))
        selected = dropdown.first_selected_option.get_attribute("value")
        log_test("Load Example sets effect measure", selected == "SMD")

        # Check threshold was set
        threshold = driver.find_element(By.ID, "clinicalThreshold").get_attribute("value")
        log_test("Load Example sets threshold", threshold == "0.5")

    except Exception as e:
        log_test("Load Example functionality", False, str(e))

def test_run_analysis(driver):
    """Test 8: Run MAFI Analysis"""
    print("\n=== TEST GROUP 8: Run Analysis ===")

    try:
        # Load example data first
        load_btn = driver.find_element(By.XPATH, "//button[contains(text(), 'Load Example')]")
        load_btn.click()
        time.sleep(0.3)

        # Click Run MAFI Analysis
        run_btn = driver.find_element(By.XPATH, "//button[contains(text(), 'Run MAFI Analysis')]")
        run_btn.click()
        time.sleep(1)

        # Check results section is visible
        results_section = driver.find_element(By.ID, "resultsSection")
        is_visible = "hidden" not in results_section.get_attribute("class")
        log_test("Results section displays", is_visible)

    except Exception as e:
        log_test("Run Analysis", False, str(e))

def test_meta_analysis_results(driver):
    """Test 9: Meta-analysis results display"""
    print("\n=== TEST GROUP 9: Meta-Analysis Results ===")

    result_elements = [
        ("pooledEffect", "Pooled Effect"),
        ("pooledCI", "95% CI"),
        ("pValue", "P-value"),
        ("numStudies", "Number of Studies"),
        ("i2Value", "I² Value"),
        ("tau2Value", "τ² Value"),
        ("qValue", "Q Statistic")
    ]

    for elem_id, name in result_elements:
        try:
            element = driver.find_element(By.ID, elem_id)
            value = element.text
            has_value = value and value != "--" and len(value) > 0
            log_test(f"{name} displayed", has_value, f"Value: {value}")
        except Exception as e:
            log_test(f"{name} displayed", False, str(e))

def test_mafi_score_display(driver):
    """Test 10: MAFI score display"""
    print("\n=== TEST GROUP 10: MAFI Score Display ===")

    try:
        mafi_value = driver.find_element(By.ID, "mafiValue").text
        has_value = mafi_value and mafi_value != "0.00"
        log_test("MAFI value displayed", has_value, f"Score: {mafi_value}")
    except Exception as e:
        log_test("MAFI value displayed", False, str(e))

    try:
        classification = driver.find_element(By.ID, "mafiClassification").text
        valid_classes = ["Robust", "Low Fragility", "Moderate Fragility", "High Fragility"]
        log_test("MAFI classification displayed", classification in valid_classes, f"Class: {classification}")
    except Exception as e:
        log_test("MAFI classification displayed", False, str(e))

    try:
        interpretation = driver.find_element(By.ID, "mafiInterpretation").text
        log_test("Interpretation displayed", len(interpretation) > 20)
    except Exception as e:
        log_test("Interpretation displayed", False, str(e))

def test_mafi_variants(driver):
    """Test 11: MAFI variants display"""
    print("\n=== TEST GROUP 11: MAFI Variants ===")

    variants = ["mafi5comp", "mafi3comp", "mafiSimple", "mafiEmpirical"]

    for variant in variants:
        try:
            score_elem = driver.find_element(By.ID, variant)
            score = score_elem.text
            has_value = score and score != "--"
            log_test(f"{variant} score displayed", has_value, f"Score: {score}")

            badge_elem = driver.find_element(By.ID, f"{variant}Badge")
            badge_text = badge_elem.text
            log_test(f"{variant} badge displayed", len(badge_text) > 0)
        except Exception as e:
            log_test(f"{variant} display", False, str(e))

def test_fragility_indices(driver):
    """Test 12: Fragility indices display"""
    print("\n=== TEST GROUP 12: Fragility Indices ===")

    indices = [
        ("dfiValue", "DFI"),
        ("sfiValue", "SFI"),
        ("cfiValue", "CFI"),
        ("effectStability", "Effect Stability")
    ]

    for elem_id, name in indices:
        try:
            element = driver.find_element(By.ID, elem_id)
            value = element.text
            log_test(f"{name} displayed", value and value != "--", f"Value: {value}")
        except Exception as e:
            log_test(f"{name} displayed", False, str(e))

    # Test progress bars
    progress_bars = ["dfiBar", "sfiBar", "cfiBar"]
    for bar_id in progress_bars:
        try:
            bar = driver.find_element(By.ID, bar_id)
            width = bar.value_of_css_property("width")
            log_test(f"{bar_id} progress bar", True)
        except Exception as e:
            log_test(f"{bar_id} progress bar", False, str(e))

def test_component_bar(driver):
    """Test 13: Component breakdown bar"""
    print("\n=== TEST GROUP 13: Component Bar ===")

    try:
        bar = driver.find_element(By.ID, "componentBar")
        segments = bar.find_elements(By.CLASS_NAME, "component-segment")
        log_test("Component bar has segments", len(segments) > 0, f"Found {len(segments)} segments")
    except Exception as e:
        log_test("Component bar", False, str(e))

    try:
        legend = driver.find_element(By.ID, "componentLegend")
        items = legend.find_elements(By.CLASS_NAME, "legend-item")
        log_test("Legend has items", len(items) >= 5, f"Found {len(items)} legend items")
    except Exception as e:
        log_test("Component legend", False, str(e))

def test_forest_plot(driver):
    """Test 14: Leave-one-out forest plot"""
    print("\n=== TEST GROUP 14: Forest Plot ===")

    try:
        forest = driver.find_element(By.ID, "forestPlot")
        rows = forest.find_elements(By.CLASS_NAME, "forest-row")
        log_test("Forest plot has rows", len(rows) > 2, f"Found {len(rows)} rows")
    except Exception as e:
        log_test("Forest plot rows", False, str(e))

    try:
        # Check for overall row
        overall_text = driver.find_element(By.XPATH, "//div[contains(@class, 'forest-row')]//strong[text()='Overall']")
        log_test("Forest plot has Overall row", overall_text is not None)
    except Exception as e:
        log_test("Forest plot Overall row", False, str(e))

    try:
        # Check for forest plot visual elements
        points = driver.find_elements(By.CLASS_NAME, "forest-point")
        log_test("Forest points displayed", len(points) > 0, f"Found {len(points)} points")

        cis = driver.find_elements(By.CLASS_NAME, "forest-ci")
        log_test("Forest CIs displayed", len(cis) > 0, f"Found {len(cis)} CIs")
    except Exception as e:
        log_test("Forest plot elements", False, str(e))

def test_grade_assessment(driver):
    """Test 15: GRADE assessment display"""
    print("\n=== TEST GROUP 15: GRADE Assessment ===")

    try:
        grade_section = driver.find_element(By.ID, "gradeAssessment")
        html_content = grade_section.get_attribute("innerHTML")
        log_test("GRADE assessment has content", len(html_content) > 100)
    except Exception as e:
        log_test("GRADE assessment", False, str(e))

    try:
        # Check for downgrade suggestion
        has_suggestion = "Suggested Downgrade" in driver.page_source or "downgrade" in driver.page_source.lower()
        log_test("GRADE downgrade suggestion present", has_suggestion)
    except Exception as e:
        log_test("GRADE downgrade suggestion", False, str(e))

def test_recommendations(driver):
    """Test 16: Recommendations display"""
    print("\n=== TEST GROUP 16: Recommendations ===")

    try:
        recs = driver.find_element(By.ID, "recommendations")
        html_content = recs.get_attribute("innerHTML")
        log_test("Recommendations section has content", len(html_content) > 100)
    except Exception as e:
        log_test("Recommendations", False, str(e))

    try:
        # Check for recommendation categories
        has_reviewers = "Systematic Reviewers" in driver.page_source
        has_developers = "Guideline Developers" in driver.page_source
        log_test("Recommendations for reviewers present", has_reviewers)
        log_test("Recommendations for guideline developers present", has_developers)
    except Exception as e:
        log_test("Recommendation categories", False, str(e))

def test_export_functions(driver):
    """Test 17: Export functionality"""
    print("\n=== TEST GROUP 17: Export Functions ===")

    try:
        export_text = driver.find_element(By.ID, "exportText")
        report = export_text.get_attribute("value")
        log_test("Report text generated", len(report) > 100, f"Report length: {len(report)}")
        log_test("Report contains MAFI", "MAFI" in report)
        log_test("Report contains leave-one-out", "Leave-One-Out" in report or "LEAVE-ONE-OUT" in report)
    except Exception as e:
        log_test("Report generation", False, str(e))

    # Test export buttons exist
    export_buttons = [
        ("Export CSV", "exportCSV"),
        ("Export JSON", "exportJSON"),
        ("Copy Report", "copyReport")
    ]

    for btn_text, func_name in export_buttons:
        try:
            btn = driver.find_element(By.XPATH, f"//button[contains(text(), '{btn_text}')]")
            log_test(f"{btn_text} button exists", btn is not None)
        except Exception as e:
            log_test(f"{btn_text} button", False, str(e))

def test_grade_tab_integration(driver):
    """Test 18: GRADE tab integration"""
    print("\n=== TEST GROUP 18: GRADE Tab Integration ===")

    try:
        # Switch to GRADE tab
        tabs = driver.find_elements(By.CLASS_NAME, "tab")
        tabs[2].click()  # GRADE tab
        time.sleep(0.3)

        # Check GRADE table exists
        grade_table = driver.find_elements(By.CLASS_NAME, "grade-table")
        log_test("GRADE table present", len(grade_table) > 0)

        # Check GRADE circles
        grade_circles = driver.find_elements(By.CLASS_NAME, "grade-circle")
        log_test("GRADE circles displayed", len(grade_circles) >= 4)

        # Check for results in GRADE tab after analysis
        grade_results = driver.find_element(By.ID, "gradeResults")
        has_results = "hidden" not in grade_results.get_attribute("class")
        log_test("GRADE results populated", has_results)

    except Exception as e:
        log_test("GRADE tab integration", False, str(e))

    # Return to calculator
    driver.find_elements(By.CLASS_NAME, "tab")[0].click()
    time.sleep(0.3)

def test_documentation_tab(driver):
    """Test 19: Documentation tab content"""
    print("\n=== TEST GROUP 19: Documentation Tab ===")

    try:
        # Switch to Documentation tab
        tabs = driver.find_elements(By.CLASS_NAME, "tab")
        tabs[3].click()  # Documentation tab
        time.sleep(0.3)

        # Check documentation content
        doc_tab = driver.find_element(By.ID, "documentationTab")
        content = doc_tab.text

        log_test("Documentation contains MAFI formula", "MAFI" in content)
        log_test("Documentation contains components", "Direction Fragility" in content)
        log_test("Documentation contains classification", "Robust" in content and "Low" in content)
        log_test("Documentation contains validation", "4,424" in content or "Cochrane" in content)
        log_test("Documentation contains references", "References" in content)

    except Exception as e:
        log_test("Documentation tab", False, str(e))

    # Return to calculator
    driver.find_elements(By.CLASS_NAME, "tab")[0].click()
    time.sleep(0.3)

def test_javascript_functions(driver):
    """Test 20: JavaScript statistical functions"""
    print("\n=== TEST GROUP 20: JavaScript Functions ===")

    # Test normalCDF
    try:
        result = driver.execute_script("return normalCDF(0)")
        log_test("normalCDF(0) = 0.5", abs(result - 0.5) < 0.01, f"Result: {result}")
    except Exception as e:
        log_test("normalCDF", False, str(e))

    try:
        result = driver.execute_script("return normalCDF(1.96)")
        log_test("normalCDF(1.96) ≈ 0.975", abs(result - 0.975) < 0.01, f"Result: {result}")
    except Exception as e:
        log_test("normalCDF(1.96)", False, str(e))

    # Test pValueFromZ
    try:
        result = driver.execute_script("return pValueFromZ(1.96)")
        log_test("pValueFromZ(1.96) ≈ 0.05", abs(result - 0.05) < 0.01, f"Result: {result}")
    except Exception as e:
        log_test("pValueFromZ", False, str(e))

    # Test gamma function
    try:
        result = driver.execute_script("return gamma(5)")
        log_test("gamma(5) = 24", abs(result - 24) < 0.1, f"Result: {result}")
    except Exception as e:
        log_test("gamma", False, str(e))

    # Test classifyMAFI
    try:
        result = driver.execute_script("return classifyMAFI(0.10)")
        log_test("classifyMAFI(0.10) = Robust", result['label'] == 'Robust', f"Result: {result}")
    except Exception as e:
        log_test("classifyMAFI robust", False, str(e))

    try:
        result = driver.execute_script("return classifyMAFI(0.60)")
        log_test("classifyMAFI(0.60) = High", result['label'] == 'High Fragility', f"Result: {result}")
    except Exception as e:
        log_test("classifyMAFI high", False, str(e))

def test_different_effect_measures(driver):
    """Test 21: Analysis with different effect measures"""
    print("\n=== TEST GROUP 21: Different Effect Measures ===")

    # Test OR (Odds Ratio)
    try:
        # Load mortality case study which uses OR
        tabs = driver.find_elements(By.CLASS_NAME, "tab")
        tabs[1].click()  # Case studies
        time.sleep(0.3)

        case_studies = driver.find_elements(By.CLASS_NAME, "case-study")
        for cs in case_studies:
            if "mortality" in (cs.get_attribute("onclick") or ""):
                cs.click()
                break
        time.sleep(0.3)

        # Run analysis
        run_btn = driver.find_element(By.XPATH, "//button[contains(text(), 'Run MAFI Analysis')]")
        run_btn.click()
        time.sleep(1)

        # Check effect measure is OR
        dropdown = Select(driver.find_element(By.ID, "effectMeasure"))
        selected = dropdown.first_selected_option.get_attribute("value")
        log_test("OR analysis runs", selected == "OR")

        # Check results are valid (OR should be around 0.7-0.9)
        pooled = driver.find_element(By.ID, "pooledEffect").text
        log_test("OR pooled effect valid", pooled != "--")

    except Exception as e:
        log_test("OR analysis", False, str(e))

def test_edge_cases(driver):
    """Test 22: Edge cases"""
    print("\n=== TEST GROUP 22: Edge Cases ===")

    # Test with minimum studies (3)
    try:
        driver.find_element(By.XPATH, "//button[contains(text(), 'Clear All')]").click()
        time.sleep(0.2)

        # Add exactly 3 studies
        rows = driver.find_elements(By.CSS_SELECTOR, "#studyTableBody tr")
        for i, row in enumerate(rows[:3]):
            inputs = row.find_elements(By.TAG_NAME, "input")
            inputs[0].clear()
            inputs[0].send_keys(f"Study {i+1}")
            inputs[1].clear()
            inputs[1].send_keys(str(-0.3 - i*0.1))
            inputs[2].clear()
            inputs[2].send_keys(str(0.1 + i*0.02))

        # Set threshold
        threshold_input = driver.find_element(By.ID, "clinicalThreshold")
        threshold_input.clear()
        threshold_input.send_keys("0.5")

        # Run analysis
        run_btn = driver.find_element(By.XPATH, "//button[contains(text(), 'Run MAFI Analysis')]")
        run_btn.click()
        time.sleep(1)

        # Check analysis completed
        results = driver.find_element(By.ID, "resultsSection")
        log_test("Analysis with 3 studies", "hidden" not in results.get_attribute("class"))

    except Exception as e:
        log_test("Minimum studies edge case", False, str(e))

    # Test insufficient studies (should show alert)
    try:
        driver.find_element(By.XPATH, "//button[contains(text(), 'Clear All')]").click()
        time.sleep(0.2)

        # Only 1 study
        row = driver.find_elements(By.CSS_SELECTOR, "#studyTableBody tr")[0]
        inputs = row.find_elements(By.TAG_NAME, "input")
        inputs[0].clear()
        inputs[0].send_keys("Single Study")
        inputs[1].clear()
        inputs[1].send_keys("-0.5")
        inputs[2].clear()
        inputs[2].send_keys("0.1")

        # Run analysis - should show alert
        run_btn = driver.find_element(By.XPATH, "//button[contains(text(), 'Run MAFI Analysis')]")
        run_btn.click()
        time.sleep(0.5)

        # Alert should appear (we just check the function doesn't crash)
        try:
            alert = driver.switch_to.alert
            alert_text = alert.text
            alert.accept()
            log_test("Insufficient studies shows alert", "3 studies" in alert_text.lower() or "least 3" in alert_text.lower())
        except:
            log_test("Insufficient studies shows alert", False, "No alert appeared")

    except Exception as e:
        log_test("Insufficient studies edge case", False, str(e))

def test_high_fragility_case(driver):
    """Test 23: High fragility case study"""
    print("\n=== TEST GROUP 23: High Fragility Case ===")

    try:
        # Load fragile case study
        tabs = driver.find_elements(By.CLASS_NAME, "tab")
        tabs[1].click()  # Case studies
        time.sleep(0.3)

        case_studies = driver.find_elements(By.CLASS_NAME, "case-study")
        for cs in case_studies:
            if "fragile" in (cs.get_attribute("onclick") or ""):
                cs.click()
                break
        time.sleep(0.3)

        # Run analysis
        run_btn = driver.find_element(By.XPATH, "//button[contains(text(), 'Run MAFI Analysis')]")
        run_btn.click()
        time.sleep(1)

        # Check classification is High or Moderate
        classification = driver.find_element(By.ID, "mafiClassification").text
        log_test("High fragility case classified correctly", "High" in classification or "Moderate" in classification, f"Class: {classification}")

        # Check interpretation box is warning/danger
        interp_box = driver.find_element(By.ID, "interpretationBox")
        box_class = interp_box.get_attribute("class")
        log_test("Interpretation shows warning/danger", "warning" in box_class or "danger" in box_class)

    except Exception as e:
        log_test("High fragility case", False, str(e))

def test_responsive_elements(driver):
    """Test 24: UI elements responsiveness"""
    print("\n=== TEST GROUP 24: UI Responsiveness ===")

    try:
        # Check grid layouts exist
        grids = driver.find_elements(By.CSS_SELECTOR, ".grid-2, .grid-3, .grid-4")
        log_test("Grid layouts present", len(grids) > 0, f"Found {len(grids)} grids")
    except Exception as e:
        log_test("Grid layouts", False, str(e))

    try:
        # Check stat boxes
        stat_boxes = driver.find_elements(By.CLASS_NAME, "stat-box")
        log_test("Stat boxes present", len(stat_boxes) >= 4, f"Found {len(stat_boxes)} stat boxes")
    except Exception as e:
        log_test("Stat boxes", False, str(e))

    try:
        # Check badges
        badges = driver.find_elements(By.CLASS_NAME, "badge")
        log_test("Badges present", len(badges) > 0, f"Found {len(badges)} badges")
    except Exception as e:
        log_test("Badges", False, str(e))

def run_all_tests():
    """Run all tests"""
    print("=" * 60)
    print("MAFI Calculator Complete - Comprehensive Selenium Test")
    print("=" * 60)

    driver = None
    try:
        driver = setup_driver()

        # Load the app
        app_path = r"C:\Users\user\OneDrive - NHS\Documents\Pairwise70\MAFI-Calculator-Complete.html"
        driver.get(f"file:///{app_path}")
        time.sleep(2)

        # Run all test groups
        test_page_load(driver)
        test_navigation_tabs(driver)
        test_study_management(driver)
        test_effect_measure_dropdown(driver)
        test_clinical_threshold_input(driver)
        test_case_studies(driver)
        test_load_example_data(driver)
        test_run_analysis(driver)
        test_meta_analysis_results(driver)
        test_mafi_score_display(driver)
        test_mafi_variants(driver)
        test_fragility_indices(driver)
        test_component_bar(driver)
        test_forest_plot(driver)
        test_grade_assessment(driver)
        test_recommendations(driver)
        test_export_functions(driver)
        test_grade_tab_integration(driver)
        test_documentation_tab(driver)
        test_javascript_functions(driver)
        test_different_effect_measures(driver)
        test_edge_cases(driver)
        test_high_fragility_case(driver)
        test_responsive_elements(driver)

    except Exception as e:
        print(f"\nCRITICAL ERROR: {e}")
        test_results["failed"] += 1
    finally:
        if driver:
            driver.quit()

    # Print summary
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    total = test_results["passed"] + test_results["failed"]
    print(f"Total Tests: {total}")
    print(f"Passed: {test_results['passed']} ({test_results['passed']/total*100:.1f}%)")
    print(f"Failed: {test_results['failed']} ({test_results['failed']/total*100:.1f}%)")

    if test_results["failed"] > 0:
        print("\nFailed Tests:")
        for test in test_results["tests"]:
            if test["status"] == "FAIL":
                print(f"  - {test['name']}: {test['details']}")

    # Save results
    results_path = r"C:\Users\user\OneDrive - NHS\Documents\Pairwise70\tests\selenium_test_results.json"
    with open(results_path, "w") as f:
        json.dump(test_results, f, indent=2)
    print(f"\nResults saved to: {results_path}")

    return test_results["failed"] == 0

if __name__ == "__main__":
    success = run_all_tests()
    exit(0 if success else 1)
