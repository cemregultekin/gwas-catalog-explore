import streamlit as st
import pandas as pd
import plotly.express as px
import requests

# --- SAYFA AYARLARI ---
st.set_page_config(page_title="GWAS Catalog Explorer", layout="wide", page_icon="🧬")

st.title("🧬 GWAS Catalog Explorer")
st.markdown("Automated lookup for the latest GWAS studies directly from OpenGWAS servers. Filter existing datasets, explore ancestry distributions, and prepare data for portability or PRS studies.")

# --- API AUTHENTICATION (SECRETS) ---
try:
    api_token = st.secrets["OPENGWAS_TOKEN"]
except Exception:
    st.error("⚠️ API configuration is missing. Please set 'OPENGWAS_TOKEN' in Streamlit Secrets.")
    st.stop()

headers = {"Authorization": f"Bearer {api_token}"}

# --- VERİ YÜKLEME (Doğrudan OpenGWAS Envanteri) ---
@st.cache_data
def load_data_from_opengwas():
    # Fetching the live inventory of OpenGWAS instead of downloading huge EBI TSV file
    res = requests.get("https://api.opengwas.io/api/gwasinfo", headers=headers)
    if res.status_code == 200:
        df = pd.DataFrame(res.json())
        
        # Standardizing column names to match your existing app logic
        df['Extract_Year'] = pd.to_numeric(df['year'], errors='coerce').fillna(2020).astype(int)
        df['N_Size'] = pd.to_numeric(df['sample_size'], errors='coerce').fillna(0).astype(int)
        return df
    return pd.DataFrame()

# Load Data
with st.spinner('Fetching active OpenGWAS core inventory...'):
    df = load_data_from_opengwas()

if df.empty:
    st.error("Failed to fetch data from OpenGWAS API. Please check your token or server status.")
    st.stop()

# --- ARAYÜZ VE MANTIK ---
if not df.empty:
    trait_col = 'trait'
    sample_col = 'population'  # Using clean population data instead of messy text string
    author_col = 'author'
    
    col1, col2 = st.columns([1, 4]) 
    
    # SOL PANEL: FİLTRELER (Tasarımın Birebir Korundu)
    with col1:
        st.header("🛠️ Filters")
        
        with st.form("filter_form"):
            st.markdown("### Disease / Trait")
            all_traits = sorted(df[trait_col].dropna().unique().tolist())
            selected_trait = st.selectbox("Select Trait", ["All"] + all_traits)
            
            st.markdown("### Ancestry")
            major_ancestries = ["European", "African", "East Asian", "South Asian", "Mixed"]
            selected_ancestries = st.multiselect("Ancestry Quick-Select", major_ancestries)
            searched_ancestry = st.text_input("Manual Ancestry Search", placeholder="e.g., Finnish, Japanese")
            
            st.markdown("### Publication Year")
            min_year_data = int(df['Extract_Year'].min()) if pd.notnull(df['Extract_Year'].min()) else 2000
            max_year_data = int(df['Extract_Year'].max()) if pd.notnull(df['Extract_Year'].max()) else 2026
            selected_year = st.slider("Minimum Year", min_year_data, max_year_data, min_year_data)
            
            st.markdown("### Data Availability")
            require_sum_stats = st.checkbox("Only show studies with Full Summary Statistics", value=True)
            
            submitted = st.form_submit_button("🚀 Apply & Visualize")
        
        if submitted:
            st.session_state["form_submitted"] = True
        
    # SAĞ PANEL: GÖRSELLEŞTİRME VE TABLOLAR (Tasarımın Birebir Korundu)
    with col2:
        st.subheader("📊 Results & Analytics")
        
        if st.session_state.get("form_submitted", False):
            # 1. Filtreleri Uygula
            res = df.copy()
            
            if selected_trait != "All":
                res = res[res[trait_col].astype(str) == selected_trait]
                
            if selected_ancestries:
                pattern = '|'.join([anc.lower() for anc in selected_ancestries])
                res = res[res[sample_col].astype(str).str.contains(pattern, case=False, na=False)]
                
            if searched_ancestry:
                res = res[res[sample_col].astype(str).str.contains(searched_ancestry, case=False, na=False)]
                
            res = res[res['Extract_Year'] >= selected_year]
            
            # 2. Sonuçları Çizdir
            if not res.empty:
                v_col1, v_col2 = st.columns(2)
                
                with v_col1:
                    # Pie Chart Logic
                    if not selected_ancestries and not searched_ancestry:
                        fig_pie = px.pie(res, names='population', values='N_Size', 
                                     title="Sample Size (N) Distribution by Ancestry",
                                     hole=0.4, color_discrete_sequence=px.colors.qualitative.Pastel)
                    else:
                        top_studies = res.nlargest(10, 'N_Size')
                        fig_pie = px.pie(top_studies, names=author_col, values='N_Size', 
                                     title="Top 10 Studies by Sample Size (N)",
                                     hole=0.4, color_discrete_sequence=px.colors.qualitative.Safe)
                        
                    st.plotly_chart(fig_pie, use_container_width=True)

                with v_col2:
                    # Line Chart Logic
                    year_counts = res['Extract_Year'].value_counts().sort_index().reset_index()
                    year_counts.columns = ['Year', 'Count']
                    fig_line = px.line(year_counts, x='Year', y='Count', 
                                       title="Studies Published Over Time",
                                       markers=True, line_shape="spline")
                    st.plotly_chart(fig_line, use_container_width=True)

                # 3. İnteraktif Tablo ve Seçim Köprüsü
                st.write(f"**Total Studies Found:** {len(res)}")
                st.info("💡 Showing first 100 rows for preview. **Select exactly TWO studies** from the table to unlock cross-population analysis.")
                
                # Selecting clean columns for the dashboard view
                display_cols = ['id', 'trait', 'population', 'sample_size', 'author', 'year']
                display_df = res[display_cols].head(100)
                
                # Multi-row selection table
                selection = st.dataframe(
                    display_df, 
                    use_container_width=True,
                    on_select="rerun",
                    selection_mode="multi-row"
                )
                
                selected_rows = selection.get("selection", {}).get("rows", [])

                if len(selected_rows) == 2:
                    study_1_id = str(display_df.iloc[selected_rows[0]]['id'])
                    study_2_id = str(display_df.iloc[selected_rows[1]]['id'])
                    
                    st.success(f"Selected Studies for Analysis: **{study_1_id}** and **{study_2_id}**")
                    
                    if st.button("🚀 Go Further Analysis & Compare", type="primary"):
                        st.session_state["auto_study_1"] = study_1_id
                        st.session_state["auto_study_2"] = study_2_id
                        st.switch_page("pages/2_⚖️_Compare_GWAS.py")
                        
                elif len(selected_rows) > 2:
                    st.warning("⚠️ Please select exactly 2 studies for comparison.")

                # Download Button
                csv_data = res[display_cols].to_csv(index=False).encode('utf-8')
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
