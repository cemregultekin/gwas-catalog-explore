import streamlit as st
import pandas as pd
import plotly.express as px
import urllib.request
import os
import re
import requests
import numpy as np

# --- PAGE CONFIGURATION ---
st.set_page_config(page_title="GWAS Catalog Explorer", layout="wide", page_icon="🧬")

st.title("🧬 GWAS Catalog Explorer")
st.markdown("Automated lookup for the latest GWAS studies directly from EBI servers. Filter existing datasets, explore ancestry distributions, and prepare data for portability or PRS studies.")

# --- SECRETS CONFIGURATION CHECK ---
try:
    api_token = st.secrets["OPENGWAS_TOKEN"]
except Exception:
    st.error("⚠️ API configuration is missing. Please set 'OPENGWAS_TOKEN' in Streamlit Secrets.")
    st.stop()

# --- DATA DOWNLOAD AND PROCESSING (Cloud Compatible) ---
@st.cache_data
def load_data():
    url = "https://www.ebi.ac.uk/gwas/api/search/downloads/studies/v1.0.2.1"
    local_file = "gwas_data.tsv"
    
    if not os.path.exists(local_file):
        with st.spinner('Downloading the latest comprehensive GWAS data from EBI servers (this happens only once)...'):
            urllib.request.urlretrieve(url, local_file)
            
    df = pd.read_csv(local_file, sep='\t', low_memory=False)
    df.columns = df.columns.str.strip()
    
    date_col = [col for col in df.columns if 'date' in col.lower()][0]
    df['Extract_Year'] = pd.to_datetime(df[date_col], errors='coerce').dt.year
    
    def get_n(text):
        nums = re.findall(r'[0-9]+(?:,[0-9]+)*', str(text))
        if nums:
            return int(nums[0].replace(',', ''))
        return 0
        
    df['N_Size'] = df['INITIAL SAMPLE SIZE'].apply(get_n)
    return df

# Load Data
try:
    df = load_data()
except Exception as e:
    st.error(f"Failed to fetch data from EBI. Error: {e}")
    st.stop()

# --- INTERFACE AND LOGIC ---
if not df.empty:
    trait_col = [col for col in df.columns if 'trait' in col.lower() or 'disease' in col.lower()][0]
    sample_col = 'INITIAL SAMPLE SIZE'
    author_col = 'FIRST AUTHOR'
    
    col1, col2 = st.columns([1, 4]) 
    
    # LEFT PANEL: FILTERS FORM
    with col1:
        st.header("🛠️ Filters")
        
        with st.form("filter_form"):
            st.markdown("### Disease / Trait")
            all_traits = sorted(df[trait_col].dropna().unique().tolist())
            selected_trait = st.selectbox("Select Trait", ["All"] + all_traits)
            
            st.markdown("### Ancestry")
            major_ancestries = ["European", "African", "East Asian", "South Asian", "Admixed American", "Hispanic", "Latino"]
            selected_ancestries = st.multiselect("Ancestry Quick-Select", major_ancestries)
            searched_ancestry = st.text_input("Manual Ancestry Search", placeholder="e.g., Finnish, Pima Indian")
            
            st.markdown("### Publication Year")
            min_year_data = int(df['Extract_Year'].min()) if pd.notnull(df['Extract_Year'].min()) else 2000
            max_year_data = int(df['Extract_Year'].max()) if pd.notnull(df['Extract_Year'].max()) else 2026
            selected_year = st.slider("Minimum Year", min_year_data, max_year_data, min_year_data)
            
            st.markdown("### Data Availability")
            require_sum_stats = st.checkbox("Only show studies with Full Summary Statistics", value=False)
            
            submitted = st.form_submit_button("🚀 Apply & Visualize")
        
        if submitted:
            st.session_state["form_submitted"] = True
        
    # RIGHT PANEL: ANALYTICS, CHARTS & GEOGRAPHY
    with col2:
        st.subheader("📊 Results & Analytics")
        
        if st.session_state.get("form_submitted", False):
            res = df.copy()
            
            if selected_trait != "All":
                res = res[res[trait_col].astype(str) == selected_trait]
                
            if selected_ancestries:
                pattern = '|'.join([anc.lower() for geom in selected_ancestries for anc in [geom]])
                res = res[res[sample_col].astype(str).str.contains(pattern, case=False, na=False)]
                
            if searched_ancestry:
                res = res[res[sample_col].astype(str).str.contains(searched_ancestry, case=False, na=False)]
                
            res = res[res['Extract_Year'] >= selected_year]
            
            if require_sum_stats:
                sum_stats_col = 'FULL SUMMARY STATISTICS'
                if sum_stats_col in res.columns:
                    res = res[res[sum_stats_col].astype(str).str.lower().str.contains('yes', na=False)]
            
            if not res.empty:
                tab_charts, tab_map = st.tabs(["📈 Statistical Charts", "🌍 Global Distribution Map"])
                
                with tab_charts:
                    v_col1, v_col2 = st.columns(2)
                    with v_col1:
                        if not selected_ancestries and not searched_ancestry:
                            res['Broad_Ancestry'] = "Other/Mixed"
                            for anc in major_ancestries:
                                res.loc[res[sample_col].str.contains(anc, case=False, na=False), 'Broad_Ancestry'] = anc
                            
                            fig_pie = px.pie(res, names='Broad_Ancestry', values='N_Size', 
                                         title="Sample Size (N) Distribution by Ancestry",
                                         hole=0.4, color_discrete_sequence=px.colors.qualitative.Pastel)
                        else:
                            top_studies = res.nlargest(10, 'N_Size')
                            fig_pie = px.pie(top_studies, names=author_col, values='N_Size', 
                                         title="Top 10 Studies by Sample Size (N)",
                                         hole=0.4, color_discrete_sequence=px.colors.qualitative.Safe)
                            
                        st.plotly_chart(fig_pie, use_container_width=True)

                    with v_col2:
                        year_counts = res['Extract_Year'].value_counts().sort_index().reset_index()
                        year_counts.columns = ['Year', 'Count']
                        fig_line = px.line(year_counts, x='Year', y='Count', 
                                           title="Studies Published Over Time",
                                           markers=True, line_shape="spline")
                        st.plotly_chart(fig_line, use_container_width=True)

                with tab_map:
                    st.markdown("##### 📍 Geographic Distribution of Study Cohorts (Global Coverage)")
                    st.caption("Dinamically resolved locations using advanced text mining on ethnic descriptions. Bubble sizes are scaled relative to Sample Size (N).")
                    
                    # 🌐 EVRENSEL ETNİSİTE VE COĞRAFYA SÖZLÜĞÜ (INITIAL SAMPLE SIZE İçin)
                    global_geo_db = {
                        # Europe & Baltics
                        'UK': {'lat': 55.3781, 'lon': -3.4360, 'lbl': 'United Kingdom'},
                        'United Kingdom': {'lat': 55.3781, 'lon': -3.4360, 'lbl': 'United Kingdom'},
                        'British': {'lat': 55.3781, 'lon': -3.4360, 'lbl': 'United Kingdom'},
                        'Finland': {'lat': 61.9241, 'lon': 25.7482, 'lbl': 'Finland'},
                        'Finnish': {'lat': 61.9241, 'lon': 25.7482, 'lbl': 'Finland'},
                        'Estonia': {'lat': 58.5953, 'lon': 25.0136, 'lbl': 'Estonia'},
                        'Estonian': {'lat': 58.5953, 'lon': 25.0136, 'lbl': 'Estonia'},
                        'Iceland': {'lat': 64.9631, 'lon': -19.0208, 'lbl': 'Iceland'},
                        'Icelandic': {'lat': 64.9631, 'lon': -19.0208, 'lbl': 'Iceland'},
                        'Sweden': {'lat': 60.1282, 'lon': 18.6435, 'lbl': 'Sweden'},
                        'Swedish': {'lat': 60.1282, 'lon': 18.6435, 'lbl': 'Sweden'},
                        'Germany': {'lat': 51.1657, 'lon': 10.4515, 'lbl': 'Germany'},
                        'German': {'lat': 51.1657, 'lon': 10.4515, 'lbl': 'Germany'},
                        'France': {'lat': 46.2276, 'lon': 2.2137, 'lbl': 'France'},
                        'French': {'lat': 46.2276, 'lon': 2.2137, 'lbl': 'France'},
                        'Netherlands': {'lat': 52.1326, 'lon': 5.2913, 'lbl': 'Netherlands'},
                        'Dutch': {'lat': 52.1326, 'lon': 5.2913, 'lbl': 'Netherlands'},
                        
                        # Asia & Middle East
                        'Korea': {'lat': 35.9078, 'lon': 127.7669, 'lbl': 'Korea'},
                        'Korean': {'lat': 35.9078, 'lon': 127.7669, 'lbl': 'Korea'},
                        'Japan': {'lat': 36.2048, 'lon': 138.2529, 'lbl': 'Japan'},
                        'Japanese': {'lat': 36.2048, 'lon': 138.2529, 'lbl': 'Japan'},
                        'China': {'lat': 35.8617, 'lon': 104.1954, 'lbl': 'China'},
                        'Chinese': {'lat': 35.8617, 'lon': 104.1954, 'lbl': 'China'},
                        'Taiwan': {'lat': 23.6978, 'lon': 120.9605, 'lbl': 'Taiwan'},
                        'Taiwanese': {'lat': 23.6978, 'lon': 120.9605, 'lbl': 'Taiwan'},
                        'Middle Eastern': {'lat': 29.2985, 'lon': 42.5510, 'lbl': 'Middle East'},
                        'Saudi': {'lat': 23.8859, 'lon': 45.0792, 'lbl': 'Saudi Arabia'},
                        'Turkey': {'lat': 38.9637, 'lon': 35.2433, 'lbl': 'Turkey'},
                        'Turkish': {'lat': 38.9637, 'lon': 35.2433, 'lbl': 'Turkey'},
                        'India': {'lat': 20.5937, 'lon': 78.9629, 'lbl': 'India'},
                        'Indian': {'lat': 20.5937, 'lon': 78.9629, 'lbl': 'India'},
                        
                        # Americas & broad fallbacks
                        'United States': {'lat': 37.0902, 'lon': -95.7129, 'lbl': 'United States'},
                        'US': {'lat': 37.0902, 'lon': -95.7129, 'lbl': 'United States'},
                        'USA': {'lat': 37.0902, 'lon': -95.7129, 'lbl': 'United States'},
                        'Hispanic': {'lat': -14.6048, 'lon': -57.6562, 'lbl': 'Hispanic/Latino'},
                        'Latino': {'lat': -14.6048, 'lon': -57.6562, 'lbl': 'Hispanic/Latino'},
                        'African': {'lat': 1.6508, 'lon': 22.5644, 'lbl': 'African (Broad)'},
                        'East Asian': {'lat': 34.0479, 'lon': 100.6197, 'lbl': 'East Asian (Broad)'},
                        'South Asian': {'lat': 22.0000, 'lon': 77.0000, 'lbl': 'South Asian (Broad)'},
                        'European': {'lat': 50.1109, 'lon': 8.6821, 'lbl': 'European (Broad)'}
                    }

                    map_rows_jittered = []
                    for idx, row in res.iterrows():
                        sample_text = str(row[sample_col])
                        n_size = row['N_Size']
                        trait = str(row[trait_col])
                        author = str(row[author_col])
                        year = int(row['Extract_Year']) if pd.notnull(row['Extract_Year']) else "Unknown"
                        accession = str(row.get('STUDY ACCESSION', 'Unknown'))
                        
                        lat, lon, resolved_name = None, None, None
                        
                        # Kelime köklerini akıllıca arayan RegEx motoru
                        for key, coord in global_geo_db.items():
                            if re.search(r'\b' + re.escape(key) + r'\b', sample_text, re.IGNORECASE):
                                lat, lon, resolved_name = coord['lat'], coord['lon'], coord['lbl']
                                break # İlk eşleşen en spesifik olanı al (Sözlük sırasına göre)

                        # Koordinat bulunduysa jittering ekle ve listeye at
                        if lat is not None:
                            jitter_amount = 1.8 
                            lat_jittered = lat + np.random.uniform(-jitter_amount, jitter_amount)
                            lon_jittered = lon + np.random.uniform(-jitter_amount, jitter_amount)
                            
                            map_rows_jittered.append({
                                'Region': resolved_name, 
                                'Latitude': lat_jittered, 
                                'Longitude': lon_jittered, 
                                'N_Size': n_size,
                                'Trait': trait,
                                'Author': author,
                                'Year': year,
                                'STUDY ACCESSION': accession
                            })

                    map_data = pd.DataFrame(map_rows_jittered)

                    if not map_data.empty:
                        map_data['size_normalized'] = np.sqrt(map_data['N_Size'].clip(lower=1)) 

                        fig_map = px.scatter_geo(
                            map_data,
                            lat="Latitude",
                            lon="Longitude",
                            size="size_normalized",
                            hover_name="Region",
                            hover_data={
                                'N_Size': ':,', 
                                'Latitude': False, 
                                'Longitude': False, 
                                'size_normalized': False,
                                'Trait': True,
                                'Author': True,
                                'Year': True,
                                'STUDY ACCESSION': True
                            },
                            size_max=16, 
                            projection="natural earth",
                            color="N_Size",
                            color_continuous_scale=px.colors.sequential.Plasma,
                            template="plotly_white",
                            opacity=0.75
                        )
                        fig_map.update_layout(margin={"r":0,"t":40,"l":0,"b":0}, geo=dict(showland=True, landcolor="#F4F4F4"))
                        st.plotly_chart(fig_map, use_container_width=True)
                    else:
                        st.info("🗺️ No geographic or ancestry-based location data could be extracted for this trait.")

                # 3. INTERACTIVE DATAFRAME & INTEGRATED PORTAL BRIDGE
                st.markdown("---")
                st.write(f"**Total Studies Found:** {len(res)}")
                st.info("💡 Showing first 100 rows for preview. **Select exactly TWO studies** from the table to unlock cross-population analysis.")
                
                display_df = res.drop(columns=['Extract_Year', 'N_Size']).head(100)
                
                selection = st.dataframe(
                    display_df, 
                    use_container_width=True,
                    on_select="rerun",
                    selection_mode="multi-row"
                )
                
                selected_rows = selection.get("selection", {}).get("rows", [])

                if len(selected_rows) == 2:
                    id_col = 'STUDY ACCESSION'
                    
                    study_1_raw = str(display_df.iloc[selected_rows[0]][id_col])
                    study_2_raw = str(display_df.iloc[selected_rows[1]][id_col])
                    
                    study_1_id = f"ebi-a-{study_1_raw}" if study_1_raw.startswith("GCST") else study_1_raw
                    study_2_id = f"ebi-a-{study_2_raw}" if study_2_raw.startswith("GCST") else study_2_raw
                    
                    st.success(f"Selected Studies: **{study_1_id}** and **{study_2_id}**")
                    
                    if st.button("🚀 Go Further Analysis & Compare", type="primary"):
                        with st.spinner("Verifying availability in OpenGWAS database..."):
                            try:
                                check_url = "https://api.opengwas.io/api/gwasinfo"
                                headers = {"Authorization": f"Bearer {api_token}", "Accept": "application/json"}
                                check_res = requests.post(check_url, json={"id": [study_1_id, study_2_id]}, headers=headers)
                                
                                if check_res.status_code == 200:
                                    available_studies = [item['id'] for item in check_res.json()]
                                    
                                    missing = []
                                    if study_1_id not in available_studies: missing.append(study_1_id)
                                    if study_2_id not in available_studies: missing.append(study_2_id)
                                    
                                    if not missing:
                                        st.session_state["auto_study_1"] = study_1_id
                                        st.session_state["auto_study_2"] = study_2_id
                                        st.switch_page("pages/2_⚖️_Compare_GWAS.py")
                                    else:
                                        st.error(f"⚠️ Unable to proceed. The following study/studies are not indexed or fully processed in OpenGWAS yet: {', '.join(missing)}. Please select different studies.")
                                else:
                                    st.error("⚠️ OpenGWAS server responded with an error. Please try again later.")
                            except Exception as e:
                                r_beta = 0 # Dummy to satisfy any inner bindings if necessary
                                st.error(f"⚠️ Connection error during verification: {e}")
                        
                elif len(selected_rows) > 2:
                    st.warning("⚠️ Please select exactly 2 studies for comparison.")

                csv_data = res.drop(columns=['Extract_Year', 'N_Size']).to_csv(index=False).encode('utf-8')
                st.download_button(
                    label="📥 Download Full Filtered Dataset (CSV)",
                    data=csv_data,
                    file_name=f"gwas_filtered_results.csv",
                    mime='text/csv',
                    type="secondary"
                )
            else:
                st.warning("⚠️ No data matches the selected criteria. Please loosen your filters.")
        else:
            st.info("👈 Please set your filters on the left and click 'Apply & Visualize' to see the results.")
