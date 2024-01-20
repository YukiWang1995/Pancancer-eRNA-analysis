import streamlit as st

st.set_page_config(
    page_title="癌症增强子数据库",
    page_icon="👋",
)

st.write("# 欢迎使用癌症增强子数据库! 👋")

st.sidebar.success("在上方选择一个演示。")

st.markdown(
    """
    癌症增强子数据库是一个专门用于计算不同癌症中增强子表达的数据库
    侧边栏是不同的应用
    ### 想了解更多吗？
    - 查看 [streamlit.io](https://streamlit.io)
    - 阅读我们的 [文档](https://docs.streamlit.io)
    - 向我们的 [邮箱](https://discuss.streamlit.io) 提问
"""
)