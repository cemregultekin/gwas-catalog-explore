import streamlit as st
import pandas as pd
import plotly.express as px
import urllib.request
import os
import re
import requests

# --- SAYFA AYARLARI ---
st.set_page_config(page_title="GWAS Catalog Explorer", layout="wide", page_icon="🧬")

st.title("🧬 GWAS Catalog Explorer")
st.markdown("Automated lookup for the latest GWAS studies directly from EBI servers. Filter existing datasets, explore ancestry distributions, and prepare data for portability or PRS studies.")

# --- SECRETS CHECK ---
try:
    api_token = st.secrets["OPENGWAS_TOKEN"]
except Exception:
    st.error("⚠️ API configuration is missing. Please set 'OPENGWAS_TOKEN' in Streamlit Secrets.")
    st.stop()

# --- VERİ İNDİRME VE İŞLEME (Bulut Uyumlu) ---
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

try:
    df = load_data()
except Exception as e:
    st.error(f"Failed to fetch data from EBI. Error: {e}")
    st.stop()

# --- ARAYÜZ VE MANTIK ---
if not df.empty:
    trait_col = [col for col in df.columns if 'trait' in col.lower() or 'disease' in col.lower()][0]
    sample_col = 'INITIAL SAMPLE SIZE'
    author_col = 'FIRST AUTHOR'
    
    col1, col2 = st.columns([1, 4]) 
    
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
                # SEKME (TABS) YAPISI EKLENDİ - Daha temiz bir görünüm için
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
                    st.markdown("##### 📍 Geographic Heatmap of Cohort Origins")
                    st.caption("Data extracted via NLP from sample descriptions. Highlights the regions where samples were collected.")
                    
                    # NLP TEXT MINING SÖZLÜĞÜ (Ülke -> ISO-3 Kod)
                    country_mapping = {
                        'UK': 'GBR', 'United Kingdom': 'GBR', 'British': 'GBR',
                        'US': 'USA', 'USA': 'USA', 'United States': 'USA', 'American': 'USA',
                        'Japan': 'JPN', 'Japanese': 'JPN', 'China': 'CHN', 'Chinese': 'CHN', 'Taiwan': 'TWN',
                        'Iceland': 'ISL', 'Icelandic': 'ISL', 'Finland': 'FIN', 'Finnish': 'FIN', 
                        'Denmark': 'DNK', 'Danish': 'DNK', 'Sweden': 'SWE', 'Swedish': 'SWE', 
                        'Norway': 'NOR', 'Norwegian': 'NOR', 'Germany': 'DEU', 'German': 'DEU', 
                        'France': 'FRA', 'French': 'FRA', 'Korea': 'KOR', 'Korean': 'KOR', 
                        'India': 'IND', 'Indian': 'IND', 'Africa': 'ZAF', 'African': 'ZAF', 
                        'Australia': 'AUS', 'Australian': 'AUS', 'Canada': 'CAN', 'Canadian': 'CAN',
                        'Netherlands': 'NLD', 'Dutch': 'NLD', 'Italy': 'ITA', 'Italian': 'ITA',
                        'Spain': 'ESP', 'Spanish': 'ESP', 'Brazil': 'BRA', 'Brazilian': 'BRA',
                        'Estonia': 'EST', 'Estonian': 'EST', 'Mexico': 'MEX', 'Mexican': 'MEX'
                    }

                    # Metinden ülke kodu çıkarma
                    res_map = res.copy()
                    res_map['ISO_Code'] = None
                    
                    for word, code in country_mapping.items():
                        # Tam kelime eşleşmesi arar (Örn: "US" bulur ama "VIRUS" içindeki US'i almaz)
                        mask = res_map[sample_col].astype(str).str.contains(r'\b' + word + r'\b', case=False, regex=True)
                        # Daha önce atanmamışsa ata
                        res_map.loc[mask & res_map['ISO_Code'].isnull(), 'ISO_Code'] = code

                    # Harita verisini grupla
                    map_data = res_map.groupby('ISO_Code')['N_Size'].sum().reset_index()

                    if not map_data.empty:
                        fig_map = px.choropleth(
                            map_data, 
                            locations="ISO_Code",
                            color="N_Size",
                            hover_name="ISO_Code",
                            color_continuous_scale=px.colors.sequential.Plasma,
                            projection="natural earth",
                        )
                        fig_map.update_layout(margin={"r":0,"t":0,"l":0,"b":0}, geo=dict(showcoastlines=True))
                        st.plotly_chart(fig_map, use_container_width=True)
                    else:
                        st.info("🗺️ No specific country data could be extracted for the selected cohort.")

                # 3. İnteraktif Tablo ve Seçim Köprüsü
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
