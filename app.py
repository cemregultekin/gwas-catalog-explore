import streamlit as st
import pandas as pd
import plotly.express as px
import urllib.request
import os
import re

# --- SAYFA AYARLARI ---
st.set_page_config(page_title="GWAS Catalog Explorer", layout="wide", page_icon="🧬")

st.title("🧬 GWAS Catalog Explorer")
st.markdown("Automated lookup for the latest GWAS studies directly from EBI servers. Filter existing datasets, explore ancestry distributions, and prepare data for portability or PRS studies.")

# --- VERİ İNDİRME VE İŞLEME (Bulut Uyumlu) ---
@st.cache_data
def load_data():
    # EBI GWAS Catalog API Güncel Linki
    url = "https://www.ebi.ac.uk/gwas/api/search/downloads/studies/v1.0.2.1"
    local_file = "gwas_data.tsv"
    
    # Sunucuda dosya yoksa internetten indirir (Sadece ilk girişte çalışır)
    if not os.path.exists(local_file):
        with st.spinner('Downloading the latest comprehensive GWAS data from EBI servers (this happens only once)...'):
            urllib.request.urlretrieve(url, local_file)
            
    # Veriyi oku
    df = pd.read_csv(local_file, sep='\t', low_memory=False)
    df.columns = df.columns.str.strip()
    
    # Yıl sütununu ayıkla
    date_col = [col for col in df.columns if 'date' in col.lower()][0]
    df['Extract_Year'] = pd.to_datetime(df[date_col], errors='coerce').dt.year
    
    # Örneklem sayısını (N) metinden saf sayıya çevir
    def get_n(text):
        nums = re.findall(r'[0-9]+(?:,[0-9]+)*', str(text))
        if nums:
            return int(nums[0].replace(',', ''))
        return 0
        
    df['N_Size'] = df['INITIAL SAMPLE SIZE'].apply(get_n)
    return df

# Veriyi Yükle
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
    
    # SOL PANEL: FİLTRELER
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
        
    # SAĞ PANEL: GÖRSELLEŞTİRME VE TABLOLAR
    with col2:
        st.subheader("📊 Results & Analytics")
        
        if submitted:
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
            
            if require_sum_stats:
                sum_stats_col = 'FULL SUMMARY STATISTICS'
                if sum_stats_col in res.columns:
                    res = res[res[sum_stats_col].astype(str).str.lower().str.contains('yes', na=False)]
            
            # 2. Sonuçları Çizdir
            if not res.empty:
                v_col1, v_col2 = st.columns(2)
                
                with v_col1:
                    # Popülasyon veya En Büyük Çalışmalar Pasta Grafiği
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
                    # Yıllara Göre Yayın Trendi Çizgi Grafiği
                    year_counts = res['Extract_Year'].value_counts().sort_index().reset_index()
                    year_counts.columns = ['Year', 'Count']
                    fig_line = px.line(year_counts, x='Year', y='Count', 
                                       title="Studies Published Over Time",
                                       markers=True, line_shape="spline")
                    st.plotly_chart(fig_line, use_container_width=True)

                # 3. İnteraktif Tablo ve İndirme
                st.write(f"**Total Studies Found:** {len(res)}")
                st.info("💡 Showing first 100 rows for preview. Use the button below to download the full dataset.")
                
                display_df = res.drop(columns=['Extract_Year', 'N_Size'])
                st.dataframe(display_df.head(100), use_container_width=True)
                
                csv_data = display_df.to_csv(index=False).encode('utf-8')
                st.download_button(
                    label="📥 Download Full Filtered Dataset (CSV)",
                    data=csv_data,
                    file_name=f"gwas_filtered_results.csv",
                    mime='text/csv',
                    type="primary"
                )
            else:
                st.warning("⚠️ No data matches the selected criteria. Please loosen your filters.")
        else:
            st.info("👈 Please set your filters on the left and click 'Apply & Visualize' to see the results.")
