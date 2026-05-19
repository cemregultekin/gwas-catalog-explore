import streamlit as st
import requests
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import time
from scipy.stats import pearsonr

st.set_page_config(page_title="Compare GWAS", page_icon="⚖️", layout="wide")

st.title("⚖️ GWAS Cross-Population Comparison")
st.markdown("Harmonize and compare Summary Statistics (Beta & EAF) between two different GWAS datasets using the OpenGWAS API.")

# --- API AUTHENTICATION (SECRETS) ---
try:
    api_token = st.secrets["OPENGWAS_TOKEN"]
except Exception:
    st.error("⚠️ API configuration is missing. Please set the 'OPENGWAS_TOKEN' in Streamlit Secrets.")
    st.stop()

# --- USER INPUT PANEL (SESSION STATE INTEGRATED) ---
if "auto_study_1" in st.session_state:
    st.session_state["id1_input"] = st.session_state["auto_study_1"]
    st.session_state["id2_input"] = st.session_state["auto_study_2"]
    del st.session_state["auto_study_1"]
    del st.session_state["auto_study_2"]

if "id1_input" not in st.session_state:
    st.session_state["id1_input"] = "ebi-a-GCST90018739"
if "id2_input" not in st.session_state:
    st.session_state["id2_input"] = "ieu-a-89"

col1, col2, col3 = st.columns(3)
with col1:
    id1 = st.text_input("Study 1 ID (Base Model)", key="id1_input")
    st.caption("Base Study (Pre-computed Top Hits will be extracted from here)")
with col2:
    id2 = st.text_input("Study 2 ID (Target Model)", key="id2_input")
    st.caption("Target Study (Full Summary Stats required for overlap lookup)")
with col3:
    snp_limit = st.number_input("Max SNPs to Compare", min_value=50, max_value=1000, value=500, step=50)

if st.button("🚀 Run Harmonization & Comparison", type="primary"):
    base_url = "https://api.opengwas.io/api"
    headers = {"Authorization": f"Bearer {api_token}", "Accept": "application/json"}

    with st.status("Running Bioinformatics Pipeline...", expanded=True) as status:
        st.write(f"1️⃣ Fetching pre-computed Top Hits from {id1}...")
        res1 = requests.post(f"{base_url}/tophits", json={"id": [id1]}, headers=headers)
        
        if res1.status_code == 200 and res1.json():
            df1 = pd.DataFrame(res1.json())
            if not df1.empty:
                # Sort by p-value and isolate to user defined limit
                df1 = df1.sort_values('p').head(snp_limit)
                snp_list = df1['rsid'].dropna().unique().tolist()
                st.write(f"✅ Isolated {len(snp_list)} significant lead SNPs. Querying targets in {id2}...")
                
                # Split rsIDs into chunks of 60 to protect payload limits
                chunks = [snp_list[i:i + 60] for i in range(0, len(snp_list), 60)]
                df2_list = []
                
                progress_bar = st.progress(0)
                for i, chunk in enumerate(chunks):
                    res2 = requests.post(f"{base_url}/associations", json={"id": [id2], "variant": chunk}, headers=headers)
                    if res2.status_code == 200 and res2.json():
                        df2_list.append(pd.DataFrame(res2.json()))
                    progress_bar.progress((i + 1) / len(chunks))
                    time.sleep(0.5)
                    
                if df2_list:
                    df2 = pd.concat(df2_list, ignore_index=True)
                    st.write(f"✅ Retrieved {len(df2)} matches from Study 2. Aligning strands...")
                    
                    # Standardize column mappings across both dataframes
                    df1 = df1.rename(columns={'beta': 'b', 'effect_allele': 'ea', 'other_allele': 'nea'})
                    df2 = df2.rename(columns={'beta': 'b', 'effect_allele': 'ea', 'other_allele': 'nea'})
                    
                    cols_to_keep = ['rsid', 'ea', 'nea', 'b', 'p', 'eaf']
                    df1_filtered = df1[[c for c in cols_to_keep if c in df1.columns]]
                    df2_filtered = df2[[c for c in cols_to_keep if c in df2.columns]]
                    
                    merged = pd.merge(df1_filtered, df2_filtered, on='rsid', suffixes=('_S1', '_S2'))

                    if not merged.empty:
                        def harmonize_all(row):
                            ea1, nea1 = str(row['ea_S1']).upper(), str(row['nea_S1']).upper()
                            b1, eaf1 = row['b_S1'], row['eaf_S1']
                            ea2, nea2 = str(row['ea_S2']).upper(), str(row['nea_S2']).upper()
                            b2, eaf2 = row['b_S2'], row['eaf_S2']
                            
                            if ea1 == ea2: 
                                return pd.Series({'harm_b_S2': b2, 'harm_eaf_S2': eaf2})
                            elif ea1 == nea2 and nea1 == ea2: 
                                return pd.Series({'harm_b_S2': -b2, 'harm_eaf_S2': 1 - eaf2})
                            else: 
                                return pd.Series({'harm_b_S2': None, 'harm_eaf_S2': None})

                        merged[['harmonized_b_S2', 'harmonized_eaf_S2']] = merged.apply(harmonize_all, axis=1)
                        merged = merged.dropna(subset=['harmonized_b_S2'])
                        
                        if not merged.empty:
                            status.update(label=f"Analysis Complete! {len(merged)} SNPs successfully harmonized.", state="complete", expanded=False)
                            
                            # --- CALCULATE STATISTICS AND RENDER CORRELATION PLOTS ---
                            r_beta, p_beta = pearsonr(merged['b_S1'], merged['harmonized_b_S2'])
                            r_eaf, p_eaf = pearsonr(merged['eaf_S1'], merged['harmonized_eaf_S2'])

                            fig = make_subplots(rows=1, cols=2, subplot_titles=(
                                f"<b>Beta Correlation</b><br><span style='font-size:12px;color:gray'>r = {r_beta:.3f} | p = {p_beta:.2e}</span>", 
                                f"<b>Allele Frequency (EAF)</b><br><span style='font-size:12px;color:gray'>r = {r_eaf:.3f} | p = {p_eaf:.2e}</span>"
                            ))

                            fig.add_trace(go.Scatter(x=merged['b_S1'], y=merged['harmonized_b_S2'], mode='markers', text=merged['rsid'], marker=dict(color='#636EFA', opacity=0.7), name='Beta'), row=1, col=1)
                            fig.add_trace(go.Scatter(x=merged['eaf_S1'], y=merged['harmonized_eaf_S2'], mode='markers', text=merged['rsid'], marker=dict(color='#EF553B', opacity=0.7), name='EAF'), row=1, col=2)

                            b_min, b_max = min(merged['b_S1'].min(), merged['harmonized_b_S2'].min()), max(merged['b_S1'].max(), merged['harmonized_b_S2'].max())
                            fig.add_shape(type="line", x0=b_min, y0=b_min, x1=b_max, y1=b_max, line=dict(color="black", dash="dash"), row=1, col=1)
                            fig.add_shape(type="line", x0=0, y0=0, x1=1, y1=1, line=dict(color="black", dash="dash"), row=1, col=2)

                            fig.update_layout(height=500, showlegend=False, template="plotly_white")
                            fig.update_xaxes(title_text="Study 1 Beta", row=1, col=1)
                            fig.update_yaxes(title_text="Study 2 Beta (Harmonized)", row=1, col=1)
                            fig.update_xaxes(title_text="Study 1 EAF", row=1, col=2)
                            fig.update_yaxes(title_text="Study 2 EAF (Harmonized)", row=1, col=2)

                            st.plotly_chart(fig, use_container_width=True)
                            
                            st.subheader("Harmonized Dataset")
                            st.dataframe(merged[['rsid', 'ea_S1', 'nea_S1', 'b_S1', 'harmonized_b_S2', 'eaf_S1', 'harmonized_eaf_S2', 'p_S1']], use_container_width=True)
                        else:
                            status.update(label="Strand flip mismatch.", state="error")
                            st.error("Overlapping SNPs were found, but they could not be reconciled due to incompatible strand/allele configurations.")
                    else:
                        status.update(label="No overlapping variants.", state="error")
                        st.error("The selected lead variants from Study 1 do not overlap with Study 2's summary statistics panel.")
                else:
                    status.update(label="No matches found.", state="error")
                    st.error("No overlapping SNPs found in Study 2. It might not contain full summary statistics.")
            else:
                status.update(label="No significant SNPs.", state="error")
                st.error("Study 1 returned an empty dataframe. No genome-wide significant top hits found.")
        else:
            status.update(label="API Error", state="error")
            st.error("Failed to fetch Top Hits from Study 1. Check ID or Token validity.")
